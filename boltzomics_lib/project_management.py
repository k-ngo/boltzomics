import os
import json
import shutil
import tempfile
import threading
from datetime import datetime
from typing import Dict, List, Optional, Union
import pandas as pd
import sys

# Directory names inside results_dir that hold app state rather than a project.
RESERVED_PROJECT_DIR_NAMES = {"_job_state"}

# Serializes read-modify-write cycles on project_metadata.json. Background queue
# workers finish concurrently, and without this two workers can each read the old
# metadata and clobber the other's results.
_METADATA_WRITE_LOCK = threading.RLock()


def _atomic_write_json(path: str, payload) -> None:
    """Write JSON via a temp file + rename so a crash can never leave a truncated
    project_metadata.json (which would make the whole project unreadable)."""
    directory = os.path.dirname(path) or "."
    os.makedirs(directory, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=directory, prefix=".tmp_", suffix=".json")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(payload, f, indent=4, default=str)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, path)
    except Exception:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except Exception:
            pass
        raise


def looks_like_project_dir(item_path: str) -> bool:
    """Heuristically decide whether a folder is a screening project.

    A project is recognized by its metadata file when present, but also by the
    artifacts a run leaves behind. Runs that crash, are cancelled, or fail every
    job never get a metadata file written, and requiring one made those projects
    invisible in the UI even though their data was on disk.
    """
    if not os.path.isdir(item_path):
        return False
    try:
        entries = os.listdir(item_path)
    except OSError:
        return False

    for entry in entries:
        if entry == "project_metadata.json":
            return True
        if entry.startswith(("screening_results_", "batch_results_")) and entry.endswith(".json"):
            return True
        if entry.startswith("boltz_results_") and os.path.isdir(os.path.join(item_path, entry)):
            return True
        if entry.startswith(("screening_", "batch_")) and entry.endswith(".yaml"):
            return True
    return False


def get_project_list(results_dir: str) -> List[str]:
    """
    Get list of existing projects from the specified results directory.

    Args:
        results_dir (str): Path to the directory containing project folders

    Returns:
        List[str]: List of project names that hold project metadata or screening artifacts
    """
    if not os.path.exists(results_dir):
        return []

    projects = []
    try:
        entries = os.listdir(results_dir)
    except OSError:
        return []

    for item in entries:
        if item in RESERVED_PROJECT_DIR_NAMES or item.startswith("."):
            continue
        item_path = os.path.join(results_dir, item)
        if looks_like_project_dir(item_path):
            projects.append(item)

    return sorted(projects)


def _collect_results_from_files(project_dir: str) -> tuple:
    """Read every screening_results_*.json in a project folder.

    Returns (results, computation_times).
    """
    all_results: List[Dict] = []
    computation_times: List[float] = []
    try:
        entries = os.listdir(project_dir)
    except OSError:
        return all_results, computation_times

    results_files = sorted(
        f for f in entries
        if f.startswith("screening_results_") and f.endswith(".json")
    )
    for results_file in results_files:
        results_path = os.path.join(project_dir, results_file)
        try:
            with open(results_path, 'r') as f:
                results_data = json.load(f)
        except Exception as exc:
            # A single corrupt shard must not hide the rest of the project.
            print(f"[WARN] Skipping unreadable results file {results_file}: {exc}")
            continue
        if isinstance(results_data, list):
            all_results.extend(results_data)
        elif isinstance(results_data, dict) and 'results' in results_data:
            all_results.extend(results_data['results'])
            if results_data.get('computation_time_seconds'):
                computation_times.append(results_data['computation_time_seconds'])
    return all_results, computation_times


def _empty_metadata(project_name: str) -> Dict:
    now = datetime.now().isoformat()
    return {
        "project_name": project_name,
        "created_date": now,
        "last_updated": now,
        "total_results": 0,
        "successful_results": 0,
        "failed_results": 0,
        "computation_time_seconds": None,
        "results": [],
    }


def ensure_project_metadata(project_name: str, results_dir: str) -> Optional[str]:
    """Create the project folder and a skeleton project_metadata.json if missing.

    Called when a screening run is queued so the project is discoverable in the UI
    from the moment work starts, rather than only after the first job succeeds.
    Returns the metadata path, or None on failure.
    """
    if not project_name:
        return None
    project_dir = os.path.join(results_dir, project_name)
    metadata_file = os.path.join(project_dir, "project_metadata.json")
    try:
        with _METADATA_WRITE_LOCK:
            os.makedirs(project_dir, exist_ok=True)
            if os.path.exists(metadata_file):
                return metadata_file
            # Reconstruct from any results shards already on disk.
            existing_results, computation_times = _collect_results_from_files(project_dir)
            metadata = _empty_metadata(project_name)
            if existing_results:
                deduped = deduplicate_results(existing_results)
                metadata["results"] = deduped
                metadata["total_results"] = len(deduped)
                metadata["successful_results"] = len(
                    [r for r in deduped if isinstance(r, dict) and r.get("status") == "Success"]
                )
                metadata["failed_results"] = metadata["total_results"] - metadata["successful_results"]
                metadata["computation_time_seconds"] = sum(computation_times) if computation_times else None
            _atomic_write_json(metadata_file, metadata)
        return metadata_file
    except Exception as exc:
        print(f"[WARN] Could not initialize project metadata for {project_name}: {exc}")
        return None

def load_project_data(project_name: str, results_dir: str) -> Optional[Dict]:
    """
    Load project data from the specified project folder, aggregating all screening_results_*.json files.
    
    Args:
        project_name (str): Name of the project to load
        results_dir (str): Path to the directory containing project folders
        
    Returns:
        Optional[Dict]: Project metadata and results, or None if loading fails
    """
    project_dir = os.path.join(results_dir, project_name)
    metadata_file = os.path.join(project_dir, "project_metadata.json")

    if not os.path.isdir(project_dir):
        return None

    try:
        metadata = None
        if os.path.exists(metadata_file):
            try:
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
            except Exception as exc:
                # Corrupt or truncated metadata: fall back to rebuilding from the
                # results shards instead of making the project unopenable.
                print(f"[WARN] project_metadata.json for {project_name} is unreadable ({exc}); rebuilding from results files.")
                metadata = None

        if metadata is None:
            # No usable metadata. The run may have been interrupted, or every job
            # may have failed before anything was committed. Rebuild what we can so
            # the project still opens and its artifacts stay reachable.
            metadata = _empty_metadata(project_name)

        # PATCH: If metadata is a list (old format), wrap in dict
        if isinstance(metadata, list):
            print("[PATCH] Detected old-format project_metadata.json (list). Wrapping in dict for compatibility.")
            metadata = {
                'project_name': project_name,
                'created_date': None,
                'last_updated': None,
                'total_results': len(metadata),
                'successful_results': len([r for r in metadata if isinstance(r, dict) and r.get("status") == "Success"]),
                'failed_results': len([r for r in metadata if isinstance(r, dict) and r.get("status") != "Success"]),
                'computation_time_seconds': None,
                'results': metadata,
                # No extra fields in old format
            }

        # Aggregate all screening_results_*.json files
        all_results, computation_times = _collect_results_from_files(project_dir)

        # Keep results the metadata already carries. Some projects (older batch
        # runs, or runs whose shards were archived) only have them here, and
        # dropping them would show an empty project.
        metadata_results = metadata.get('results') if isinstance(metadata, dict) else None
        if isinstance(metadata_results, list) and metadata_results:
            all_results = metadata_results + all_results

        # Deduplicate results
        deduped_results = deduplicate_results(all_results)
        metadata['results'] = deduped_results
        # Keep the summary counters consistent with what we actually loaded, so a
        # rebuilt project does not report stale totals.
        metadata['total_results'] = len(deduped_results)
        metadata['successful_results'] = len(
            [r for r in deduped_results if isinstance(r, dict) and r.get("status") == "Success"]
        )
        metadata['failed_results'] = metadata['total_results'] - metadata['successful_results']
        # Optionally, sum computation times if desired
        metadata['computation_time_seconds'] = sum(computation_times) if computation_times else None
        # --- Load extra fields if present ---
        for extra_field in ["template_cif_path", "binding_pocket_constraints", "boltz_commands"]:
            if extra_field in metadata:
                pass  # already present
            else:
                metadata[extra_field] = None
        return metadata
    except Exception as e:
        print(f"Error loading project data: {str(e)}")
        return None

def delete_project(project_name: str, results_dir: str) -> bool:
    """
    Delete a project folder and all its contents.
    
    Args:
        project_name (str): Name of the project to delete
        results_dir (str): Path to the directory containing project folders
        
    Returns:
        bool: True if deletion was successful, False otherwise
    """
    project_dir = os.path.join(results_dir, project_name)
    
    if not os.path.exists(project_dir):
        return False
    
    try:
        shutil.rmtree(project_dir)
        return True
    except Exception as e:
        print(f"Error deleting project: {str(e)}")
        return False

def deduplicate_results(results: List[Dict]) -> List[Dict]:
    """
    Remove duplicate entries from results, keeping the most complete and recent entries.
    
    Args:
        results: List of prediction result dictionaries
        
    Returns:
        List of deduplicated results
    """
    if not results:
        return results

    # Drop anything that is not a dict; a malformed entry in one shard must not
    # take down the whole project load.
    results = [r for r in results if isinstance(r, dict)]
    if not results:
        return []

    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(results)

    # Create a composite key for identifying duplicates. Entries written by the
    # force-update/recovery path use 'protein'/'drug', so accept those as aliases
    # and tolerate rows missing the fields entirely.
    def _key_series(primary: str, alias: str):
        if primary in df.columns:
            series = df[primary]
            if alias in df.columns:
                series = series.fillna(df[alias])
        elif alias in df.columns:
            series = df[alias]
        else:
            series = pd.Series([""] * len(df), index=df.index)
        return series.fillna("").astype(str)

    df['composite_key'] = _key_series('protein_name', 'protein') + '|' + _key_series('drug_name', 'drug')

    # Group by composite key and keep the best entry
    deduplicated_results = []
    
    for key, group in df.groupby('composite_key'):
        if len(group) == 1:
            # No duplicates, keep as is
            result = group.iloc[0].to_dict()
            result.pop('composite_key', None)  # Remove the temporary key
            deduplicated_results.append(result)
        else:
            # Multiple entries for same protein-drug combination
            # Score each entry based on completeness and recency
            best_entry = None
            best_score = -1
            
            for _, row in group.iterrows():
                score = 0
                
                # Completeness score: more non-null values = higher score
                numeric_fields = ['ic50_um', 'pic50', 'affinity_probability', 'confidence', 'ptm', 'iptm', 'avg_plddt']
                for field in numeric_fields:
                    if field in row and pd.notna(row[field]):
                        score += 1
                
                # Status score: Success > other statuses
                if row.get('status') == 'Success':
                    score += 10
                
                # Recency score: prefer entries with workspace/design info (indicating newer runs)
                if row.get('workspace') and row.get('design'):
                    score += 5
                
                # If scores are equal, prefer the first one (original order)
                if score > best_score:
                    best_score = score
                    best_entry = row
            
            if best_entry is not None:
                result = best_entry.to_dict()
                result.pop('composite_key', None)  # Remove the temporary key
                deduplicated_results.append(result)
    
    return deduplicated_results

def save_screening_results(results: Union[List[Dict], Dict], 
                      project_name: str, 
                      results_dir: str,
                      computation_time: Optional[float] = None,
                      template_cif_path: Optional[str] = None,
                      binding_pocket_constraints: Optional[dict] = None,
                      boltz_commands: Optional[list] = None) -> Optional[str]:
    """
    Save screening results to project-specific folder.
    
    Args:
        results (Union[List[Dict], Dict]): Results to save, either as list or dict with metadata
        project_name (str): Name of the project
        results_dir (str): Path to the directory containing project folders
        computation_time (Optional[float]): Computation time in seconds
        
    Returns:
        Optional[str]: Path to saved results file, or None if saving failed
    """
    # Background queue workers call this concurrently; the metadata update below is
    # a read-modify-write cycle and must not interleave.
    _METADATA_WRITE_LOCK.acquire()
    try:
        # Create project directory
        project_dir = os.path.join(results_dir, project_name)
        if not os.path.exists(project_dir):
            os.makedirs(project_dir)

        # Create filename with timestamp. Seconds resolution is not enough when
        # several workers finish together, so keep adding a suffix until the name
        # is free instead of silently overwriting another worker's shard.
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        filename = f"screening_results_{timestamp}.json"
        filepath = os.path.join(project_dir, filename)
        suffix = 0
        while os.path.exists(filepath):
            suffix += 1
            filename = f"screening_results_{timestamp}_{suffix}.json"
            filepath = os.path.join(project_dir, filename)

        # Extract results list
        results_list = results if isinstance(results, list) else results.get('results', [])

        # Add computation time to results metadata
        results_with_metadata = {
            "computation_time_seconds": computation_time,
            "timestamp": timestamp,
            "results": results_list
        }

        # Save results
        _atomic_write_json(filepath, results_with_metadata)

        # Load existing project metadata if it exists
        metadata_file = os.path.join(project_dir, "project_metadata.json")
        existing_metadata = {}
        if os.path.exists(metadata_file):
            try:
                with open(metadata_file, 'r') as f:
                    existing_metadata = json.load(f)
            except:
                existing_metadata = {}
        
        # Get existing results if available
        existing_results = []
        if 'results' in existing_metadata:
            existing_results = existing_metadata['results']
        elif os.path.exists(metadata_file):
            # Try to load from the most recent results file
            results_files = [f for f in os.listdir(project_dir) if f.startswith("screening_results_") and f.endswith(".json")]
            if results_files:
                results_files.sort(reverse=True)
                latest_results_file = os.path.join(project_dir, results_files[0])
                try:
                    with open(latest_results_file, 'r') as f:
                        results_data = json.load(f)
                        # Handle both old and new formats
                        if isinstance(results_data, list):
                            existing_results = results_data
                        elif isinstance(results_data, dict) and 'results' in results_data:
                            existing_results = results_data['results']
                        else:
                            existing_results = []
                except:
                    existing_results = []
        
        # Combine existing and new results
        all_results = existing_results + results_list
        
        # Deduplicate the combined results
        deduplicated_results = deduplicate_results(all_results)
        
        # Update project metadata
        metadata = {
            "project_name": project_name,
            "created_date": existing_metadata.get("created_date", datetime.now().isoformat()),
            "last_updated": datetime.now().isoformat(),
            "total_results": len(deduplicated_results),
            "successful_results": len([r for r in deduplicated_results if isinstance(r, dict) and r.get("status") == "Success"]),
            "failed_results": len([r for r in deduplicated_results if not isinstance(r, dict) or r.get("status") != "Success"]),
            "computation_time_seconds": computation_time,
            "results": deduplicated_results,
            "template_cif_path": template_cif_path if template_cif_path is not None else existing_metadata.get("template_cif_path"),
            "binding_pocket_constraints": binding_pocket_constraints if binding_pocket_constraints is not None else existing_metadata.get("binding_pocket_constraints"),
            "boltz_commands": boltz_commands if boltz_commands is not None else existing_metadata.get("boltz_commands")
        }
        
        # Handle backward compatibility - preserve existing computation time if not provided
        if computation_time is None and 'computation_time_seconds' in existing_metadata:
            metadata["computation_time_seconds"] = existing_metadata["computation_time_seconds"]
        elif computation_time is None:
            # Set to None if no computation time available (backward compatibility)
            metadata["computation_time_seconds"] = None
        
        # Save updated metadata atomically so an interrupted write cannot leave a
        # truncated file that makes the project unreadable on the next page load.
        _atomic_write_json(metadata_file, metadata)

        return filepath
    except Exception as e:
        print(f"Error saving results: {str(e)}")
        return None
    finally:
        _METADATA_WRITE_LOCK.release()

def rename_results_in_project(project_name: str, old_name: str, new_name: str, rename_type: str, results_dir: str) -> bool:
    """
    Rename protein or drug names in all project files and associated folders.
    
    Args:
        project_name (str): Name of the project
        old_name (str): Old name to replace
        new_name (str): New name to use
        rename_type (str): Either "protein" or "drug"
        results_dir (str): Path to the directory containing project folders
        
    Returns:
        bool: True if renaming was successful, False otherwise
    """
    try:
        project_dir = os.path.join(results_dir, project_name)
        if not os.path.exists(project_dir):
            print(f"Project directory not found: {project_dir}")
            return False
        
        # Validate rename_type
        if rename_type not in ["protein", "drug"]:
            print(f"Invalid rename_type: {rename_type}. Must be 'protein' or 'drug'")
            return False
        
        # Get the field name to rename
        field_name = "protein_name" if rename_type == "protein" else "drug_name"
        
        # Track if any changes were made
        changes_made = False
        
        # 1. Update all screening_results_*.json files
        results_files = [f for f in os.listdir(project_dir) if f.startswith("screening_results_") and f.endswith(".json")]
        for results_file in results_files:
            results_path = os.path.join(project_dir, results_file)
            try:
                with open(results_path, 'r') as f:
                    results_data = json.load(f)
                
                # Handle both old and new formats
                results_list = results_data if isinstance(results_data, list) else results_data.get('results', [])
                
                # Check if any results need renaming
                file_changed = False
                for result in results_list:
                    if isinstance(result, dict) and result.get(field_name) == old_name:
                        result[field_name] = new_name
                        file_changed = True
                
                if file_changed:
                    # Update the results in the data structure
                    if isinstance(results_data, list):
                        results_data = results_list
                    else:
                        results_data['results'] = results_list
                    
                    # Write back to file
                    with open(results_path, 'w') as f:
                        json.dump(results_data, f, indent=4, default=str)
                    
                    changes_made = True
                    print(f"Updated {results_file}")
                    
            except Exception as e:
                print(f"Error updating {results_file}: {str(e)}")
                continue
        
        # 2. Update project_metadata.json
        metadata_file = os.path.join(project_dir, "project_metadata.json")
        if os.path.exists(metadata_file):
            try:
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                # Handle both old and new metadata formats
                if isinstance(metadata, dict) and 'results' in metadata:
                    metadata_changed = False
                    for result in metadata['results']:
                        if isinstance(result, dict) and result.get(field_name) == old_name:
                            result[field_name] = new_name
                            metadata_changed = True
                    
                    if metadata_changed:
                        with open(metadata_file, 'w') as f:
                            json.dump(metadata, f, indent=4, default=str)
                        changes_made = True
                        print(f"Updated project_metadata.json")
                elif isinstance(metadata, list):
                    # Old format metadata
                    metadata_changed = False
                    for result in metadata:
                        if isinstance(result, dict) and result.get(field_name) == old_name:
                            result[field_name] = new_name
                            metadata_changed = True
                    
                    if metadata_changed:
                        with open(metadata_file, 'w') as f:
                            json.dump(metadata, f, indent=4, default=str)
                        changes_made = True
                        print(f"Updated project_metadata.json")
                        
            except Exception as e:
                print(f"Error updating project_metadata.json: {str(e)}")
        
        # 3. Rename associated Boltz result folders
        # Look for folders that contain the old name in their design name
        boltz_folders = [f for f in os.listdir(project_dir) if f.startswith("boltz_results_")]
        for folder in boltz_folders:
            try:
                # Extract the design name from the folder name
                # Format: boltz_results_<workspace>_<design>
                if folder.startswith("boltz_results_"):
                    parts = folder.split("_", 2)  # Split into ["boltz", "results", "workspace_design"]
                    if len(parts) >= 3:
                        workspace_design = parts[2]
                        # Check if the design name contains the old name
                        if old_name in workspace_design:
                            # Create new design name by replacing old_name with new_name
                            new_design_name = workspace_design.replace(old_name, new_name)
                            new_folder_name = f"boltz_results_{parts[1]}_{new_design_name}"
                            
                            old_folder_path = os.path.join(project_dir, folder)
                            new_folder_path = os.path.join(project_dir, new_folder_name)
                            
                            # Rename the folder
                            if os.path.exists(old_folder_path) and not os.path.exists(new_folder_path):
                                os.rename(old_folder_path, new_folder_path)
                                changes_made = True
                                print(f"Renamed folder: {folder} -> {new_folder_name}")
                            
                            # Also rename the YAML file if it exists
                            yaml_file = f"{parts[1]}_{workspace_design}.yaml"
                            new_yaml_file = f"{parts[1]}_{new_design_name}.yaml"
                            yaml_path = os.path.join(project_dir, yaml_file)
                            new_yaml_path = os.path.join(project_dir, new_yaml_file)
                            
                            if os.path.exists(yaml_path) and not os.path.exists(new_yaml_path):
                                os.rename(yaml_path, new_yaml_path)
                                changes_made = True
                                print(f"Renamed YAML file: {yaml_file} -> {new_yaml_file}")
                                
            except Exception as e:
                print(f"Error renaming folder {folder}: {str(e)}")
                continue
        
        if changes_made:
            print(f"Successfully renamed '{old_name}' to '{new_name}' in project '{project_name}'")
            return True
        else:
            print(f"No instances of '{old_name}' found in project '{project_name}'")
            return True  # Still return True as this is not an error
            
    except Exception as e:
        print(f"Error during rename operation: {str(e)}")
        return False

def get_active_workspace_data():
    """Returns the dictionary for the currently active workspace."""
    import streamlit as st
    active_id = st.session_state.active_workspace_id
    return st.session_state.workspaces.get(active_id, {})

def save_workspace_to_file(workspace_name, workspace_data):
    """Saves a specific workspace to its JSON file."""
    import os
    import json
    from datetime import datetime
    import numpy as np
    import streamlit as st
    def sanitize_filename(name):
        import re
        safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
        safe_name = re.sub(r'_+', '_', safe_name)
        safe_name = safe_name.strip('_')
        if not safe_name:
            safe_name = "workspace"
        return safe_name
    def get_workspace_file_path(workspace_name):
        script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        safe_name = sanitize_filename(workspace_name)
        return os.path.join(script_dir, f"workspace_{safe_name}.json")
    filepath = get_workspace_file_path(workspace_name)
    if not filepath or filepath.strip() == "":
        print("Error: Invalid file path for workspace save")
        st.error("Failed to save workspace data: Invalid file path")
        st.session_state['last_save_success'] = False
        return
    try:
        data_to_save = {
            "name": workspace_data.get("name", workspace_name),
            "ligand_smiles": workspace_data.get("ligand_smiles", ""),
            "protein_sequence": workspace_data.get("protein_sequence", ""),
            "prediction_history": workspace_data.get("prediction_history", []),
            "current_prediction": workspace_data.get("current_prediction"),
            "last_saved": datetime.now().isoformat()
        }
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                if hasattr(obj, 'dtype') and obj.dtype.names:
                    result = []
                    for row in obj:
                        row_dict = {}
                        for name in obj.dtype.names:
                            value = row[name]
                            if isinstance(value, (np.integer, np.floating, np.bool_)):
                                row_dict[name] = value.item()
                            else:
                                row_dict[name] = value
                        result.append(row_dict)
                    return result
                else:
                    return obj.tolist()
            elif isinstance(obj, (np.integer, np.floating, np.bool_)):
                return obj.item()
            elif isinstance(obj, dict):
                return {key: convert_numpy(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(item) for item in obj]
            else:
                return obj
        data_to_save = convert_numpy(data_to_save)
        directory = os.path.dirname(filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
        temp_filepath = filepath + '.tmp'
        with open(temp_filepath, 'w') as f:
            json.dump(data_to_save, f, indent=4)
        if os.path.exists(filepath):
            os.replace(temp_filepath, filepath)
        else:
            os.rename(temp_filepath, filepath)
        st.session_state['last_save_success'] = True
    except Exception as e:
        st.error(f"Failed to save workspace data: {str(e)}")
        st.session_state['last_save_success'] = False
        print(f"Error saving workspace to {filepath}: {e}")
