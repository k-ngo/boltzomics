"""
Improved Mutation Discovery Module
Simplified version inspired by CATVariant approach but streamlined for drug screening workflow
"""

import requests
import re
from typing import Dict, List, Optional, Tuple
import streamlit as st
import pandas as pd

# Amino acid property groups for categorizing mutations
AA_PROPERTIES = {
    'hydrophobic': {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'},
    'polar': {'S', 'T', 'N', 'Q', 'Y', 'C'},
    'positively_charged': {'K', 'R', 'H'},
    'negatively_charged': {'D', 'E'},
    'special': {'G', 'P'}  # Glycine and Proline have special conformational properties
}

# BLOSUM62-based severity scoring for common substitutions
BLOSUM_SEVERITY = {
    # Conservative substitutions (low severity)
    ('I', 'V'): 0.1, ('V', 'I'): 0.1, ('L', 'I'): 0.2, ('I', 'L'): 0.2,
    ('M', 'L'): 0.2, ('L', 'M'): 0.2, ('F', 'Y'): 0.2, ('Y', 'F'): 0.2,
    ('K', 'R'): 0.1, ('R', 'K'): 0.1, ('D', 'E'): 0.1, ('E', 'D'): 0.1,
    ('S', 'T'): 0.1, ('T', 'S'): 0.1, ('N', 'Q'): 0.1, ('Q', 'N'): 0.1,

    # Semi-conservative substitutions (moderate severity)
    ('A', 'G'): 0.3, ('G', 'A'): 0.3, ('A', 'S'): 0.3, ('S', 'A'): 0.3,
    ('L', 'F'): 0.4, ('F', 'L'): 0.4, ('I', 'M'): 0.4, ('M', 'I'): 0.4,
    ('K', 'Q'): 0.5, ('Q', 'K'): 0.5, ('R', 'Q'): 0.5, ('Q', 'R'): 0.5,

    # Radical substitutions (high severity)
    ('G', 'P'): 0.8, ('P', 'G'): 0.8, ('P', 'R'): 0.8, ('R', 'P'): 0.8,
    ('C', 'F'): 0.9, ('F', 'C'): 0.9, ('W', 'C'): 0.9, ('C', 'W'): 0.9,
}

def find_closest_protein_match(sequence: str) -> Dict:
    """
    Identify protein using multi-strategy approach:
    1. Check if input is a gene name/accession (e.g., "CYP2B6", "P20813")
    2. Length-based search for exact/near-exact matches (fast)
    3. BLAST search for fragments/subsequences (comprehensive)
    """
    try:
        # Check if input looks like a gene name or accession (short, no sequence characters)
        trimmed = sequence.strip()
        if len(trimmed) < 50 and ' ' not in trimmed:
            # Try as gene name or accession
            keyword_result = search_by_keyword(trimmed)
            if keyword_result and "error" not in keyword_result:
                return keyword_result

        # Clean sequence
        clean_seq = re.sub(r'[^A-Z]', '', sequence.upper())
        if len(clean_seq) < 10:
            return {"error": "Sequence too short for reliable identification (minimum 10 amino acids)"}

        seq_length = len(clean_seq)
        url = "https://rest.uniprot.org/uniprotkb/search"

        # Strategy 1: Try exact length match first (fastest for full sequences)
        query = f'reviewed:true AND length:[{seq_length} TO {seq_length}]'
        params = {
            'query': query,
            'format': 'json',
            'size': 50,
            'fields': 'accession,protein_name,organism_name,length,sequence,gene_names'
        }

        response = requests.get(url, params=params, timeout=20)
        if response.status_code == 200:
            data = response.json()
            results = data.get('results', [])

            # Look for exact sequence match
            for result in results:
                seq = result.get('sequence', {}).get('value', '')
                if seq == clean_seq:
                    return {
                        "accession": result.get('primaryAccession', ''),
                        "protein_name": extract_protein_name(result),
                        "organism": result.get('organism', {}).get('scientificName', ''),
                        "length": len(seq),
                        "similarity": 1.0,
                        "sequence": seq,
                        "gene_name": extract_gene_name(result),
                        "search_strategy": "Exact Match"
                    }

        # Strategy 2: Try with length tolerance
        # For potential fragments, search LONGER proteins (user seq might be a fragment)
        # For full sequences, also check slightly shorter/longer variants
        tolerance = min(50, seq_length // 10)
        query = f'reviewed:true AND length:[{max(1, seq_length-tolerance)} TO {seq_length+200}]'  # Search up to +200aa to find full proteins
        params = {
            'query': query,
            'format': 'json',
            'size': 100,
            'fields': 'accession,protein_name,organism_name,length,sequence,gene_names'
        }

        response = requests.get(url, params=params, timeout=20)
        if response.status_code == 200:
            data = response.json()
            results = data.get('results', [])

            best_match = None
            best_similarity = 0

            for result in results:
                seq = result.get('sequence', {}).get('value', '')
                if not seq:
                    continue

                # Calculate simple sequence identity
                if len(seq) == seq_length:
                    # Same length - direct comparison
                    similarity = sum(1 for a, b in zip(clean_seq, seq) if a == b) / seq_length
                else:
                    # Different lengths - check if one is contained in the other
                    if clean_seq in seq:
                        # User sequence is a perfect subsequence
                        similarity = 0.95
                    elif seq in clean_seq:
                        # Database sequence is contained in user sequence
                        similarity = 0.90
                    else:
                        # Partial match on overlapping region
                        min_len = min(len(clean_seq), len(seq))
                        similarity = sum(1 for a, b in zip(clean_seq[:min_len], seq[:min_len]) if a == b) / min_len

                if similarity > best_similarity and similarity >= 0.8:  # Require 80% identity minimum
                    best_similarity = similarity
                    best_match = result

            if best_match:
                return {
                    "accession": best_match.get('primaryAccession', ''),
                    "protein_name": extract_protein_name(best_match),
                    "organism": best_match.get('organism', {}).get('scientificName', ''),
                    "length": best_match.get('sequence', {}).get('length', 0),
                    "similarity": best_similarity,
                    "sequence": best_match.get('sequence', {}).get('value', ''),
                    "gene_name": extract_gene_name(best_match),
                    "search_strategy": "Length-based Search"
                }

        # Strategy 3: BLAST search (fallback for fragments/subsequences)
        # This is the key for finding full proteins when given fragments
        blast_result = blast_search_uniprot(clean_seq)
        if blast_result and "error" not in blast_result:
            return blast_result

        return {"error": f"No protein match found with ≥80% sequence identity"}

    except Exception as e:
        return {"error": f"Protein search failed: {str(e)}"}

def blast_search_uniprot(sequence: str) -> Dict:
    """
    Use NCBI BLAST to find proteins, then map to UniProt.
    This handles sequences that exist in NCBI but may not be in UniProt's reviewed entries.

    Strategy:
    1. BLAST against NCBI nr database
    2. Extract protein identifiers (gene names, protein names)
    3. Search UniProt using those identifiers
    """
    try:
        import time
        import re

        # Use EBI's NCBI BLAST service
        ebi_blast_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"

        blast_params = {
            'email': 'noreply@example.com',
            'program': 'blastp',
            'stype': 'protein',
            'sequence': sequence,
            'database': 'uniprotkb_swissprot',  # Use Swiss-Prot for curated proteins
            'alignments': 10,
            'scores': 10
        }

        # Submit job
        response = requests.post(ebi_blast_url, data=blast_params, timeout=30)
        if response.status_code != 200:
            return {"error": f"BLAST submission failed: {response.status_code}"}

        job_id = response.text.strip()

        # Poll for results
        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
        max_wait = 90  # Allow more time for NCBI BLAST
        start_time = time.time()

        while time.time() - start_time < max_wait:
            status_response = requests.get(status_url, timeout=10)
            if status_response.status_code == 200:
                status = status_response.text.strip()

                if status == 'FINISHED':
                    # Get results
                    results_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/out"
                    results_response = requests.get(results_url, timeout=30)

                    if results_response.status_code == 200:
                        # Parse BLAST output and map to UniProt
                        return parse_blast_and_map_to_uniprot(results_response.text, sequence)
                    else:
                        return {"error": "Failed to retrieve BLAST results"}

                elif status == 'FAILURE' or status == 'ERROR':
                    return {"error": "BLAST search failed"}

                elif status == 'RUNNING':
                    time.sleep(4)
                else:
                    time.sleep(3)
            else:
                return {"error": f"BLAST status check failed: {status_response.status_code}"}

        return {"error": "BLAST search timed out"}

    except Exception as e:
        return {"error": f"BLAST search error: {str(e)}"}

def parse_blast_and_map_to_uniprot(blast_output: str, query_seq: str) -> Dict:
    """
    Parse BLAST output to extract UniProt accessions and return best match.

    Strategy:
    1. Extract UniProt accessions directly from BLAST hits (SP:P12345 format)
    2. Fetch full protein data from UniProt for the top hit
    3. Return best match without arbitrary similarity thresholds
    """
    try:
        import re

        # Pattern to extract UniProt accessions from BLAST output
        # Format: "SP:P20813" or just "P20813"
        uniprot_accession_pattern = r'SP:([A-Z0-9]{6,10})|^>([A-Z0-9]{6,10})'

        # Find all UniProt accessions in the BLAST output
        accessions = []
        for line in blast_output.split('\n'):
            # Look for lines with "SP:" prefix (Swiss-Prot entries)
            if 'SP:' in line:
                match = re.search(r'SP:([A-Z0-9]{6,10})', line)
                if match:
                    accessions.append(match.group(1))

        # If no accessions found, try alternative parsing
        if not accessions:
            # Look for accession-like patterns in the first significant hit
            matches = re.findall(r'\b([A-Z][0-9][A-Z0-9]{3,4}[0-9])\b', blast_output[:3000])
            accessions = [m for m in matches if len(m) == 6]

        if not accessions:
            return {"error": "No UniProt accessions found in BLAST results"}

        # Try each accession, starting with the first (best hit)
        for accession in accessions[:5]:  # Check top 5 hits
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

            try:
                response = requests.get(uniprot_url, timeout=10)
                if response.status_code == 200:
                    result = response.json()

                    # Extract sequence
                    seq = result.get('sequence', {}).get('value', '')
                    if not seq:
                        continue

                    # Calculate actual similarity for reporting
                    query_len = len(query_seq)
                    seq_len = len(seq)

                    # Find best alignment position (handles fragments and partial matches)
                    if query_len == seq_len:
                        # Same length - direct comparison
                        similarity = sum(1 for a, b in zip(query_seq, seq) if a == b) / query_len
                    else:
                        # Different lengths - find best sliding window match
                        best_matches = 0
                        best_alignment_len = 0

                        if query_len < seq_len:
                            # Query is shorter - slide it along the database sequence
                            for i in range(seq_len - query_len + 1):
                                matches = sum(1 for a, b in zip(query_seq, seq[i:i+query_len]) if a == b)
                                if matches > best_matches:
                                    best_matches = matches
                                    best_alignment_len = query_len
                        else:
                            # Database sequence is shorter - slide it along query
                            for i in range(query_len - seq_len + 1):
                                matches = sum(1 for a, b in zip(query_seq[i:i+seq_len], seq) if a == b)
                                if matches > best_matches:
                                    best_matches = matches
                                    best_alignment_len = seq_len

                        similarity = best_matches / best_alignment_len if best_alignment_len > 0 else 0.0

                    # Return the best BLAST hit regardless of similarity threshold
                    return {
                        "accession": result.get('primaryAccession', accession),
                        "protein_name": extract_protein_name(result),
                        "organism": result.get('organism', {}).get('scientificName', ''),
                        "length": seq_len,
                        "similarity": similarity,
                        "sequence": seq,
                        "gene_name": extract_gene_name(result),
                        "search_strategy": "BLAST Search"
                    }
            except Exception as e:
                continue  # Try next accession

        return {"error": "Could not retrieve protein data from UniProt"}

    except Exception as e:
        return {"error": f"BLAST parsing error: {str(e)}"}

def search_by_keyword(keyword: str) -> Dict:
    """
    Search UniProt by gene name, accession, or protein name.
    Examples: "CYP2B6", "P20813", "cytochrome P450 2B6"
    """
    try:
        url = "https://rest.uniprot.org/uniprotkb/search"

        # Try as accession first
        if re.match(r'^[A-Z0-9]{6,10}$', keyword.upper()):
            params = {
                'query': f'accession:{keyword.upper()}',
                'format': 'json',
                'size': 1,
                'fields': 'accession,protein_name,organism_name,length,sequence,gene_names'
            }
            response = requests.get(url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                results = data.get('results', [])
                if results:
                    result = results[0]
                    return {
                        "accession": result.get('primaryAccession', ''),
                        "protein_name": extract_protein_name(result),
                        "organism": result.get('organism', {}).get('scientificName', ''),
                        "length": result.get('sequence', {}).get('length', 0),
                        "similarity": 1.0,
                        "sequence": result.get('sequence', {}).get('value', ''),
                        "gene_name": extract_gene_name(result),
                        "search_strategy": "Accession Lookup"
                    }

        # Try as gene name
        params = {
            'query': f'gene:{keyword} AND reviewed:true',
            'format': 'json',
            'size': 5,
            'fields': 'accession,protein_name,organism_name,length,sequence,gene_names'
        }
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            results = data.get('results', [])
            if results:
                # Prefer human proteins
                for result in results:
                    org = result.get('organism', {}).get('scientificName', '')
                    if 'sapiens' in org.lower():
                        return {
                            "accession": result.get('primaryAccession', ''),
                            "protein_name": extract_protein_name(result),
                            "organism": org,
                            "length": result.get('sequence', {}).get('length', 0),
                            "similarity": 1.0,
                            "sequence": result.get('sequence', {}).get('value', ''),
                            "gene_name": extract_gene_name(result),
                            "search_strategy": "Gene Name Lookup"
                        }
                # Return first result if no human found
                result = results[0]
                return {
                    "accession": result.get('primaryAccession', ''),
                    "protein_name": extract_protein_name(result),
                    "organism": result.get('organism', {}).get('scientificName', ''),
                    "length": result.get('sequence', {}).get('length', 0),
                    "similarity": 1.0,
                    "sequence": result.get('sequence', {}).get('value', ''),
                    "gene_name": extract_gene_name(result),
                    "search_strategy": "Gene Name Lookup"
                }

        return {"error": f"No matches found for keyword: {keyword}"}

    except Exception as e:
        return {"error": f"Keyword search failed: {str(e)}"}

def calculate_enhanced_similarity(seq1: str, seq2: str) -> float:
    """Enhanced sequence similarity calculation with multiple scoring factors.
    Returns a value between 0.0 and 1.0."""
    if not seq1 or not seq2:
        return 0.0

    len1, len2 = len(seq1), len(seq2)
    min_len, max_len = min(len1, len2), max(len1, len2)

    if min_len < 10:
        return 0.0

    # For very different lengths, check if shorter sequence is contained in longer one
    if len1 != len2:
        shorter, longer = (seq1, seq2) if len1 < len2 else (seq2, seq1)
        if shorter in longer:
            # Perfect subsequence match - high but not perfect score
            # Position matters, so central matches score higher
            pos = longer.index(shorter)
            position_score = 1.0 - (abs(pos - (len(longer) - len(shorter))/2) / len(longer))
            return 0.85 + (position_score * 0.10)  # 0.85-0.95 range

    # Calculate identity score with alignment consideration
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    identity = matches / min_len

    # Length similarity factor
    length_similarity = 1.0 - (abs(len1 - len2) / max_len)

    # Conservative vs radical substitutions weighting
    conservative_score = 0
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            if is_conservative_substitution(seq1[i], seq2[i]):
                conservative_score += 0.8
            else:
                conservative_score += 0.2
        else:
            conservative_score += 1.0

    # FIX: Divide by min_len (the actual number of positions scored), not max_len
    conservative_similarity = conservative_score / min_len

    # Combined score - weight identity heavily, then conservative substitutions
    # Length difference is less important for finding the right protein
    final_similarity = (identity * 0.6 + conservative_similarity * 0.3 + length_similarity * 0.1)

    # Ensure score is between 0 and 1
    return min(1.0, max(0.0, final_similarity))

def is_conservative_substitution(res1: str, res2: str) -> bool:
    """Check if amino acid substitution is conservative."""
    # Define conservative groups
    conservative_groups = [
        {'I', 'L', 'V', 'M', 'A', 'G'},  # Hydrophobic/aliphatic
        {'F', 'Y', 'W'},                  # Aromatic
        {'K', 'R', 'H'},                  # Positively charged
        {'D', 'E'},                       # Negatively charged
        {'S', 'T', 'N', 'Q'},             # Polar uncharged
        {'C'},                            # Special (sulfur)
        {'P'},                            # Special (proline)
    ]

    for group in conservative_groups:
        if res1 in group and res2 in group:
            return True
    return False

def search_by_similarity(sequence: str, seq_length: int) -> List[Dict]:
    """Fallback similarity-based search using fragment matching.
    This helps find full proteins when the user provides a subsequence."""
    try:
        # Create query fragments - use larger fragments for better specificity
        fragment_length = min(25, max(15, seq_length // 3))
        fragments = [sequence[i:i+fragment_length] for i in range(0, len(sequence)-fragment_length+1, fragment_length//2)]

        url = "https://rest.uniprot.org/uniprotkb/search"
        similarity_results = []

        # Use multiple fragments to find proteins that contain the user sequence
        for fragment in fragments[:5]:  # Use more fragments for better coverage
            params = {
                'query': f'sequence:"{fragment}" AND reviewed:true',
                'format': 'json',
                'size': 20,  # Get more results to increase chance of finding the right protein
                'fields': 'accession,protein_name,organism_name,length,sequence,gene_names,go'
            }

            try:
                response = requests.get(url, params=params, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    results = data.get('results', [])

                    for result in results:
                        seq = result.get('sequence', {}).get('value', '')
                        if seq and len(seq) >= seq_length:  # Only consider proteins >= user sequence length
                            # Check if user sequence is contained in this protein
                            if sequence in seq:
                                # Perfect subsequence match
                                similarity_score = 0.95
                            else:
                                # Calculate similarity for the overlapping region
                                similarity_score = calculate_enhanced_similarity(sequence, seq)

                            if similarity_score > 0.3:  # Keep matches with decent similarity
                                result_info = {
                                    'result': result,
                                    'strategy': 'Fragment Search',
                                    'similarity': similarity_score,
                                    'accession': result.get('primaryAccession', ''),
                                    'protein_name': extract_protein_name(result),
                                    'organism': result.get('organism', {}).get('scientificName', '').lower(),
                                    'gene_name': extract_gene_name(result)
                                }
                                similarity_results.append(result_info)

            except Exception:
                continue

        # Remove duplicates and sort by similarity
        unique_results = {}
        for result_info in similarity_results:
            acc = result_info['accession']
            if acc not in unique_results or result_info['similarity'] > unique_results[acc]['similarity']:
                unique_results[acc] = result_info

        return sorted(unique_results.values(), key=lambda x: x['similarity'], reverse=True)[:10]

    except Exception:
        return []

def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """Calculate sequence similarity using local alignment approach."""
    if not seq1 or not seq2:
        return 0.0

    # For different lengths, focus on the overlapping region
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))

    if min_len < 10:
        return 0.0

    # Calculate identity in the overlapping region
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    identity = matches / min_len

    # Apply length penalty
    length_penalty = abs(len(seq1) - len(seq2)) / max_len
    similarity = identity * (1 - length_penalty * 0.5)

    return similarity

def extract_protein_name(result: Dict) -> str:
    """Extract protein name from UniProt result."""
    name_info = result.get('proteinDescription', {})
    rec_name = name_info.get('recommendedName', {}).get('fullName', {}).get('value', '')

    if rec_name:
        return rec_name

    # Fallback to other names
    alt_names = name_info.get('alternativeNames', [])
    if alt_names:
        return alt_names[0].get('fullName', {}).get('value', 'Unknown')

    return "Unknown Protein"

def extract_gene_name(result: Dict) -> str:
    """Extract gene name from UniProt result."""
    genes = result.get('genes', [])
    if genes:
        preferred_genes = [g.get('geneName', {}).get('value', '') for g in genes if g.get('geneName', {}).get('value')]
        return preferred_genes[0] if preferred_genes else ""
    return ""

def get_protein_variants(accession: str) -> List[Dict]:
    """
    Enhanced variant retrieval using CATVariant-inspired multi-database approach.
    Integrates EBI Proteins, UniProt, and cross-reference databases for comprehensive variant collection.
    Also enriches variants with SIFT, PolyPhen, and REVEL prediction scores where available.
    """
    all_variants = []

    # Database 1: EBI Proteins Variation API (primary source like CATVariant)
    ebi_variants = get_ebi_variants(accession)
    all_variants.extend(ebi_variants)

    # Database 2: UniProt variant features (fallback/backup)
    uniprot_variants = get_uniprot_variants(accession)
    all_variants.extend(uniprot_variants)

    # Remove duplicates based on position and change
    unique_variants = deduplicate_variants(all_variants)

    # Enrich variants with prediction scores (SIFT, PolyPhen, REVEL) from EBI
    unique_variants = enrich_variants_with_predictions(unique_variants, accession)

    # Sort by prediction scores if available, then evidence count
    unique_variants.sort(key=lambda x: (
        x.get('has_predictions', False),
        x.get('revel_score', 0),
        x.get('evidence_count', 0)
    ), reverse=True)

    return unique_variants  # Return all found variants (unlimited)

def get_ebi_variants(accession: str) -> List[Dict]:
    """
    Retrieve variants from EBI Proteins Variation API with rich annotations.
    This mirrors CATVariant's primary variant source.
    """
    try:
        url = f"https://www.ebi.ac.uk/proteins/api/variation/{accession}"

        headers = {'Accept': 'application/json'}
        response = requests.get(url, headers=headers, timeout=20)

        if response.status_code != 200:
            return []

        data = response.json()
        features = data.get('features', [])

        variants = []
        for feature in features:
            if feature.get('type') == 'VARIANT':
                variant_info = parse_ebi_variant_feature(feature)
                if variant_info:
                    variant_info['source'] = 'EBI Proteins'
                    variants.append(variant_info)

        return variants

    except Exception as e:
        # Silently fall back to UniProt
        return []

def get_uniprot_variants(accession: str) -> List[Dict]:
    """
    Retrieve variants from UniProt as fallback database.
    """
    try:
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'accession:{accession}',
            'format': 'json',
            'size': 1,
            'fields': 'ft_variant,sequence'
        }

        response = requests.get(url, params=params, timeout=15)
        if response.status_code != 200:
            return []

        data = response.json()
        results = data.get('results', [])

        if not results:
            return []

        result = results[0]
        sequence = result.get('sequence', {}).get('value', '')
        features = result.get('features', [])

        variants = []
        for feature in features:
            if feature.get('type') in ['Natural variant', 'Sequence variant', 'Mutagenesis']:
                variant_info = parse_variant_feature(feature, sequence)
                if variant_info:
                    variant_info['source'] = 'UniProt'
                    variants.append(variant_info)

        return variants

    except Exception as e:
        print(f"UniProt API error: {e}")
        return []

def parse_ebi_variant_feature(feature: Dict) -> Optional[Dict]:
    """
    Parse EBI variant feature with rich clinical annotations and prediction scores.
    Extracts SIFT, PolyPhen, and other pathogenicity predictions if available.
    """
    try:
        begin = feature.get('begin', '')
        end = feature.get('end', '')
        # EBI API uses 'wildType' and 'mutatedType' fields
        wild_type = feature.get('wildType', feature.get('wild_type', ''))
        # Try both 'mutatedType', 'alternativeSequence', and 'alternative_aa'
        alt_aa = feature.get('mutatedType', feature.get('alternativeSequence', feature.get('alternative_aa', '')))

        if not begin or not wild_type or not alt_aa:
            return None

        # Handle multiple alternative amino acids
        alt_aas = alt_aa.split(',') if ',' in alt_aa else [alt_aa]

        variants = []
        for alt in alt_aas:
            if alt and alt != '?':
                # Extract clinical significance from association
                associations = feature.get('association', [])
                clinical_significance = [assoc.get('name', '') for assoc in associations if assoc.get('disease')]

                # Get consequence type
                consequence = feature.get('consequenceType', '')

                # Extract prediction scores from EBI data
                predictions = feature.get('predictions', [])
                sift_score = None
                sift_prediction = None
                polyphen_score = None
                polyphen_prediction = None

                for pred in predictions:
                    pred_type = pred.get('predAlgorithmNameType', '').lower()
                    if 'sift' in pred_type:
                        sift_score = pred.get('score')
                        sift_prediction = pred.get('predictionValType', '')
                    elif 'polyphen' in pred_type:
                        polyphen_score = pred.get('score')
                        polyphen_prediction = pred.get('predictionValType', '')

                variant_info = {
                    'position': int(begin),
                    'wt_residue': wild_type,
                    'mut_residue': alt,
                    'description': f"{wild_type}{begin}{alt}",
                    'clinical_significance': ', '.join(clinical_significance) if clinical_significance else '',
                    'consequence': consequence,
                    'evidence_count': len(feature.get('xrefs', [])),  # Use xrefs count as evidence
                    'variant_type': 'Natural variant',
                    'source_db': 'EBI',
                    'xrefs': feature.get('xrefs', []),
                    # Prediction scores
                    'sift_score': sift_score,
                    'sift_prediction': sift_prediction,
                    'polyphen_score': polyphen_score,
                    'polyphen_prediction': polyphen_prediction,
                    'has_predictions': sift_score is not None or polyphen_score is not None
                }
                variants.append(variant_info)

        return variants[0] if variants else None

    except Exception:
        return None

def extract_cross_references(result: Dict, db_name: str) -> List[str]:
    """Extract cross-references for specific database."""
    try:
        xrefs = result.get('references', [])
        return [ref.get('id', '') for ref in xrefs if db_name.lower() in ref.get('source', '').lower()]
    except Exception:
        return []

def enrich_variants_with_predictions(variants: List[Dict], accession: str) -> List[Dict]:
    """
    Enrich variants with prediction scores from additional sources if not already present.
    This is a lightweight enrichment that only adds scores for variants that lack them.
    """
    # For now, return variants as-is since EBI already provides SIFT/PolyPhen
    # In the future, could add calls to:
    # - dbNSFP for REVEL, CADD, MutationAssessor, etc.
    # - Ensembl VEP API for comprehensive predictions
    # - AlphaMissense for AI-based predictions

    for variant in variants:
        if not variant.get('has_predictions'):
            # Use our heuristic-based severity if no predictions available
            wt_res = variant.get('wt_residue', '')
            mut_res = variant.get('mut_residue', '')

            if wt_res and mut_res:
                change_type = categorize_amino_acid_change(wt_res, mut_res)
                severity = predict_severity(wt_res, mut_res, change_type)
                variant['heuristic_severity'] = severity
                variant['heuristic_change_type'] = change_type

    return variants

def deduplicate_variants(variants: List[Dict]) -> List[Dict]:
    """
    Remove duplicate variants based on position and amino acid change.
    Keep the one with highest evidence count.
    """
    seen_variants = {}
    unique_variants = []

    for variant in variants:
        key = (variant.get('position'), variant.get('wt_residue'), variant.get('mut_residue'))

        if key not in seen_variants:
            seen_variants[key] = variant
            unique_variants.append(variant)
        else:
            # Keep the variant with more evidence or richer annotations
            existing = seen_variants[key]
            if (variant.get('evidence_count', 0) > existing.get('evidence_count', 0) or
                variant.get('source') == 'EBI Proteins'):  # Prefer EBI data
                unique_variants[unique_variants.index(existing)] = variant
                seen_variants[key] = variant

    return unique_variants

def parse_variant_feature(feature: Dict, sequence: str) -> Optional[Dict]:
    """Parse a variant feature from UniProt data using the alternativeSequence field."""
    try:
        location = feature.get('location', {})
        start = location.get('start', {}).get('value')
        description = feature.get('description', [''])[0] if feature.get('description') else ''
        evidence = feature.get('evidences', [])
        alt_sequence = feature.get('alternativeSequence', '')

        if not start:
            return None

        # Extract mutation information from alternativeSequence
        wt_residue = None
        mut_residue = None

        # Handle alternativeSequence structure
        if isinstance(alt_sequence, dict):
            wt_residue = alt_sequence.get('originalSequence', '')
            alt_sequences = alt_sequence.get('alternativeSequences', [])
            mut_residue = alt_sequences[0] if alt_sequences else ''
        elif isinstance(alt_sequence, str) and len(alt_sequence) == 1:
            # Simple case: just the mutated residue
            mut_residue = alt_sequence
            # Get wild-type from sequence
            if start <= len(sequence):
                wt_residue = sequence[start - 1]
        elif alt_sequence:
            # Try to parse from string
            mut_residue = alt_sequence
            if start <= len(sequence):
                wt_residue = sequence[start - 1]

        # Validate position is within sequence
        if start <= len(sequence) and wt_residue and mut_residue:
            actual_wt = sequence[start - 1] if start > 0 else ''

            # Extract dbSNP references if available
            dbsnp_refs = []
            for ref in feature.get('featureCrossReferences', []):
                ref_id = ref.get('id', '')
                if 'rs' in ref_id:
                    dbsnp_refs.append(ref_id)

            return {
                'position': int(start),
                'wt_residue': wt_residue,
                'mut_residue': mut_residue,
                'description': description or f"{wt_residue}{start}{mut_residue}",
                'evidence_count': len(evidence),
                'variant_type': feature.get('type', 'Variant'),
                'actual_wt_residue': actual_wt,
                'dbsnp_refs': dbsnp_refs
            }

    except Exception as e:
        print(f"Error parsing variant: {e}")

    return None

def categorize_amino_acid_change(wt_res: str, mut_res: str) -> str:
    """Categorize the type of amino acid change."""
    if wt_res == mut_res:
        return "Synonymous"

    # Find property groups for each residue
    wt_groups = []
    mut_groups = []

    for group, residues in AA_PROPERTIES.items():
        if wt_res in residues:
            wt_groups.append(group)
        if mut_res in residues:
            mut_groups.append(group)

    # Check for changes involving special residues
    if wt_res in {'G', 'P'} or mut_res in {'G', 'P'}:
        if wt_res in {'G', 'P'} and mut_res in {'G', 'P'}:
            return "Special-to-Special"
        else:
            return "Special-to-Regular"

    # Check for charge changes
    wt_charged = any(group in wt_groups for group in ['positively_charged', 'negatively_charged'])
    mut_charged = any(group in mut_groups for group in ['positively_charged', 'negatively_charged'])

    if wt_charged and mut_charged:
        return "Charge-to-Charge"
    elif wt_charged or mut_charged:
        return "Charge-to-Neutral"

    # Check for hydrophobic to polar changes
    wt_hydro = 'hydrophobic' in wt_groups
    mut_hydro = 'hydrophobic' in mut_groups
    wt_polar = 'polar' in wt_groups
    mut_polar = 'polar' in mut_groups

    if wt_hydro and mut_polar:
        return "Hydrophobic-to-Polar"
    elif wt_polar and mut_hydro:
        return "Polar-to-Hydrophobic"
    elif wt_hydro and mut_hydro:
        return "Hydrophobic-to-Hydrophobic"
    elif wt_polar and mut_polar:
        return "Polar-to-Polar"

    return "Other"

def predict_severity(wt_res: str, mut_res: str, change_type: str) -> str:
    """Predict functional severity of amino acid change."""
    if wt_res == mut_res:
        return "None"

    # Use BLOSUM-based scoring if available
    blosum_score = BLOSUM_SEVERITY.get((wt_res, mut_res), 0.5)

    # Adjust based on change type
    if "Synonymous" in change_type:
        return "None"
    elif "Special" in change_type:
        return "High"
    elif "Charge" in change_type:
        if blosum_score < 0.2:
            return "Low"
        elif blosum_score < 0.5:
            return "Moderate"
        else:
            return "High"
    elif "Hydrophobic-to-Polar" in change_type or "Polar-to-Hydrophobic" in change_type:
        return "Moderate"
    elif blosum_score < 0.3:
        return "Low"
    elif blosum_score < 0.6:
        return "Moderate"
    else:
        return "High"

def predict_structural_impact(wt_res: str, mut_res: str, change_type: str, severity: str) -> str:
    """Predict structural impact based on amino acid properties."""
    if wt_res == mut_res:
        return "No structural change"

    # Special cases for structurally important residues
    if wt_res == 'P' or mut_res == 'P':
        return "Likely affects backbone conformation (proline)"
    elif wt_res == 'G' or mut_res == 'G':
        return "Likely affects flexibility (glycine)"
    elif wt_res in {'C', 'W'} or mut_res in {'C', 'W'}:
        return "May affect disulfide bonds or packing"

    # Based on severity
    if severity == "High":
        return "Likely disrupts secondary structure or core packing"
    elif severity == "Moderate":
        return "May affect local structure or surface properties"
    else:
        return "Minimal structural impact expected"

def assess_mutation_druggability_impact(wt_res: str, mut_res: str, position: int, protein_length: int) -> Dict:
    """Assess how mutation might impact drug binding."""
    # Simple heuristic based on position and residue properties
    relative_position = position / protein_length

    impact_assessment = {
        'binding_site_likelihood': 'Unknown',
        'pocket_change_potential': 'Unknown',
        'druggability_impact': 'Unknown'
    }

    # N-terminal and C-terminal regions are less likely to be in binding pockets
    if relative_position < 0.1 or relative_position > 0.9:
        impact_assessment['binding_site_likelihood'] = 'Low'
    elif 0.3 <= relative_position <= 0.7:
        impact_assessment['binding_site_likelihood'] = 'Medium-High'

    # Check for changes that might affect binding pockets
    if wt_res in AA_PROPERTIES['hydrophobic'] and mut_res in AA_PROPERTIES['polar']:
        impact_assessment['pocket_change_potential'] = 'High (hydrophobic to polar)'
    elif wt_res in AA_PROPERTIES['polar'] and mut_res in AA_PROPERTIES['hydrophobic']:
        impact_assessment['pocket_change_potential'] = 'High (polar to hydrophobic)'
    elif 'Charge' in categorize_amino_acid_change(wt_res, mut_res):
        impact_assessment['pocket_change_potential'] = 'High (charge change)'

    # Overall druggability impact
    if impact_assessment['binding_site_likelihood'] == 'Medium-High' and 'High' in impact_assessment['pocket_change_potential']:
        impact_assessment['druggability_impact'] = 'High - likely affects drug binding'
    elif impact_assessment['pocket_change_potential'] != 'Unknown':
        impact_assessment['druggability_impact'] = 'Moderate - may affect drug binding'
    else:
        impact_assessment['druggability_impact'] = 'Low - unlikely to affect drug binding'

    return impact_assessment

def discover_mutations_for_sequence(sequence: str, max_mutations: int = 20, chain_start: int = 1) -> List[Dict]:
    """
    Enhanced mutation discovery that maps database variants to user-provided sequence region.
    Even if the user sequence is a subset of the full wild-type protein, only mutations
    within the user sequence region will be displayed.
    """
    # Step 1: Find closest protein match
    protein_info = find_closest_protein_match(sequence)

    if "error" in protein_info:
        st.error(f"Protein identification failed: {protein_info['error']}")
        return []

    # Display protein match information
    similarity = protein_info.get('similarity', 0)
    user_seq_length = len(re.sub(r'[^A-Z]', '', sequence.upper()))
    full_protein_length = protein_info.get('length', 0)

    if similarity >= 0.9:
        st.success(f"Found exact match: {protein_info['protein_name']} ({protein_info['accession']})")
    elif similarity >= 0.7:
        st.info(f"Found close match: {protein_info['protein_name']} ({protein_info['accession']}) - {similarity:.1%} similarity")
    else:
        st.warning(f"Found partial match: {protein_info['protein_name']} ({protein_info['accession']})")

    # Show sequence mapping information
    if user_seq_length < full_protein_length:
        st.info(f"**Sequence Mapping**: Your sequence ({user_seq_length} aa) is shorter than the full protein ({full_protein_length} aa). Only mutations within your sequence region will be shown.")

    # Step 2: Get variants for this protein
    variants = get_protein_variants(protein_info['accession'])

    if not variants:
        st.warning("No known variants found for this protein.")
        return []

    # Step 3: Map variants to user sequence region
    processed_mutations = []
    clean_user_sequence = re.sub(r'[^A-Z]', '', sequence.upper())
    clean_full_sequence = protein_info.get('sequence', '')

    # Find mapping between user sequence and full protein
    mapping_info = map_user_sequence_to_full_protein(clean_user_sequence, clean_full_sequence)

    if mapping_info['start_position'] == -1:
        st.error("Could not map your sequence to the reference protein. Please check sequence format.")
        return []

    user_start = mapping_info['start_position']
    user_end = mapping_info['end_position']
    coverage = mapping_info['coverage']
    identity = mapping_info.get('identity', coverage)  # Fallback to coverage for backward compatibility

    st.info(f"**Mapping**: Your sequence corresponds to positions {user_start}-{user_end} of the full protein ({coverage:.1%} coverage)")

    # Process all variants and filter to user sequence region
    for variant in variants:  # Process all available variants
        canonical_pos = variant['position']
        wt_res = variant['wt_residue']
        mut_res = variant['mut_residue']

        # Only include variants within user sequence region
        if user_start <= canonical_pos <= user_end:
            # Calculate position in user sequence (1-based)
            position_in_user = canonical_pos - user_start + 1

            # Verify wild-type residue matches at this position in the user sequence
            if position_in_user <= len(clean_user_sequence):
                actual_wt = clean_user_sequence[position_in_user - 1] if position_in_user > 0 else ''

                # Allow some mismatches if the mapping isn't perfect, but prioritize exact matches
                # Get the similarity score from protein_info
                similarity_score = protein_info.get('similarity', 0)
                if actual_wt.upper() == wt_res.upper() or similarity_score < 0.8:
                    # Prefer database prediction scores over heuristics
                    if variant.get('has_predictions'):
                        # Use SIFT/PolyPhen predictions from database
                        sift_pred = variant.get('sift_prediction', '')
                        polyphen_pred = variant.get('polyphen_prediction', '')

                        # Map predictions to severity
                        if 'deleterious' in sift_pred.lower() or 'damaging' in polyphen_pred.lower():
                            severity = 'High'
                        elif 'tolerated' in sift_pred.lower() or 'benign' in polyphen_pred.lower():
                            severity = 'Low'
                        else:
                            severity = 'Moderate'

                        change_type = f"SIFT: {sift_pred}, PolyPhen: {polyphen_pred}" if sift_pred or polyphen_pred else "Predicted"
                    else:
                        # Fall back to heuristic predictions
                        change_type = variant.get('heuristic_change_type') or categorize_amino_acid_change(wt_res, mut_res)
                        severity = variant.get('heuristic_severity') or predict_severity(wt_res, mut_res, change_type)

                    structural_impact = predict_structural_impact(wt_res, mut_res, change_type, severity)
                    druggability_impact = assess_mutation_druggability_impact(wt_res, mut_res, position_in_user, len(clean_user_sequence))

                    # User-friendly numbering starts from chain_start
                    display_position = position_in_user + chain_start - 1

                    mutation_data = {
                        'mutation_code': f"{wt_res}{display_position}{mut_res}",
                        'canonical_code': f"{wt_res}{canonical_pos}{mut_res}",
                        'position_in_input': display_position,
                        'canonical_position': canonical_pos,
                        'position_in_user_sequence': position_in_user,
                        'user_sequence_start': user_start,
                        'user_sequence_end': user_end,
                        'wt_residue': wt_res,
                        'mut_residue': mut_res,
                        'change_type': change_type,
                        'severity': severity,
                        'structural_impact': structural_impact,
                        'druggability_impact': druggability_impact['druggability_impact'],
                        'evidence_count': variant['evidence_count'],
                        'variant_type': variant['variant_type'],
                        'description': variant['description'],
                        'uniprot_accession': protein_info['accession'],
                        'protein_name': protein_info['protein_name'],
                        'similarity': similarity,
                        'sequence_match': actual_wt.upper() == wt_res.upper(),
                        # Include prediction scores
                        'sift_score': variant.get('sift_score'),
                        'sift_prediction': variant.get('sift_prediction'),
                        'polyphen_score': variant.get('polyphen_score'),
                        'polyphen_prediction': variant.get('polyphen_prediction'),
                        'has_predictions': variant.get('has_predictions', False),
                        # Include literature references
                        'xrefs': variant.get('xrefs', []),
                        'clinical_significance': variant.get('clinical_significance', '')
                    }

                    processed_mutations.append(mutation_data)

    # Sort by severity (High > Moderate > Low) and evidence count
    severity_order = {'High': 3, 'Moderate': 2, 'Low': 1, 'None': 0}
    processed_mutations.sort(key=lambda x: (severity_order.get(x['severity'], 0), x['evidence_count']), reverse=True)

    # Return all found mutations (unlimited)
    return processed_mutations

def map_user_sequence_to_full_protein(user_seq: str, full_seq: str) -> Dict:
    """
    Map user sequence to full protein sequence to find the corresponding region.
    Returns mapping information including start/end positions and coverage.
    """
    if not user_seq or not full_seq:
        return {'start_position': -1, 'end_position': -1, 'coverage': 0.0}

    user_len = len(user_seq)
    full_len = len(full_seq)

    # If sequences are the same length, do direct comparison
    if user_len == full_len:
        if user_seq == full_seq:
            return {'start_position': 1, 'end_position': full_len, 'coverage': 1.0}
        else:
            # Try to find best alignment
            best_match = find_best_alignment(user_seq, full_seq)
            return best_match

    # If user sequence is shorter, find best matching region in full protein
    best_match = {'start_position': -1, 'end_position': -1, 'coverage': 0.0, 'identity': 0.0}

    for start_pos in range(0, full_len - user_len + 1):
        segment = full_seq[start_pos:start_pos + user_len]
        identity = sum(1 for i in range(user_len) if user_seq[i] == segment[i])
        identity_percent = identity / user_len

        if identity_percent > best_match['identity']:
            # Coverage = how much of the full protein is covered by user sequence
            coverage_percent = user_len / full_len
            best_match = {
                'start_position': start_pos + 1,  # Convert to 1-based
                'end_position': start_pos + user_len,
                'coverage': coverage_percent,  # % of full protein covered
                'identity': identity_percent    # % of matching amino acids
            }

        # Early exit if we find a perfect match
        if identity_percent >= 0.95:
            break

    return best_match

def find_best_alignment(seq1: str, seq2: str) -> Dict:
    """Find the best alignment between two sequences of similar length."""
    len_diff = abs(len(seq1) - len(seq2))

    # If length difference is small, do simple offset alignment
    if len_diff <= 10:
        if len(seq1) <= len(seq2):
            # seq1 is shorter or equal
            for offset in range(len_diff + 1):
                segment = seq2[offset:offset + len(seq1)]
                identity = sum(1 for i in range(len(seq1)) if seq1[i] == segment[i])
                identity_percent = identity / len(seq1)
                if identity_percent >= 0.7:
                    coverage_percent = len(seq1) / len(seq2)
                    return {
                        'start_position': offset + 1,
                        'end_position': offset + len(seq1),
                        'coverage': coverage_percent,
                        'identity': identity_percent
                    }
        else:
            # seq2 is shorter
            for offset in range(len_diff + 1):
                segment = seq1[offset:offset + len(seq2)]
                identity = sum(1 for i in range(len(seq2)) if seq2[i] == segment[i])
                identity_percent = identity / len(seq2)
                if identity_percent >= 0.7:
                    coverage_percent = len(seq2) / len(seq1)
                    return {
                        'start_position': offset + 1,
                        'end_position': offset + len(seq2),
                        'coverage': coverage_percent,
                        'identity': identity_percent
                    }

    return {'start_position': -1, 'end_position': -1, 'coverage': 0.0, 'identity': 0.0}

def format_mutations_for_display(mutations: List[Dict]) -> pd.DataFrame:
    """Format mutations for display in Streamlit with mapping information."""
    if not mutations:
        return pd.DataFrame()

    df_data = []
    for mut in mutations:
        # Create detailed mutation description
        canonical_pos = mut.get('canonical_position', mut.get('position_in_input', 'Unknown'))
        user_pos = mut.get('position_in_user_sequence', 'Unknown')

        if user_pos != 'Unknown' and canonical_pos != 'Unknown':
            position_info = f"User:{user_pos} → Ref:{canonical_pos}"
        else:
            position_info = f"Pos {canonical_pos}"

        df_data.append({
            "Mutation": mut['mutation_code'],
            "Reference": mut['canonical_code'],
            "Position Mapping": position_info,
            "Change Type": mut['change_type'],
            "Severity": mut['severity'],
            "Structural Impact": mut['structural_impact'],
            "Drug Binding Impact": mut['druggability_impact'],
            "Evidence Count": mut['evidence_count'],
            "Protein": mut['protein_name'][:30] + "..." if len(mut['protein_name']) > 30 else mut['protein_name']
        })

    return pd.DataFrame(df_data)

def format_mutations_for_input(mutations: List[Dict], max_count: int = 10) -> str:
    """Format mutations for input to mutation screening."""
    if not mutations:
        return ""

    # Take top mutations by severity and evidence
    selected_mutations = mutations[:max_count]
    mutation_codes = [m['mutation_code'] for m in selected_mutations]

    return ", ".join(mutation_codes)

def display_mutation_summary_statistics(mutations: List[Dict]):
    """Display summary statistics for discovered mutations."""
    if not mutations:
        return

    # Count by severity
    severity_counts = {}
    change_type_counts = {}

    for mut in mutations:
        severity = mut.get('severity', 'Unknown')
        change_type = mut.get('change_type', 'Unknown')

        severity_counts[severity] = severity_counts.get(severity, 0) + 1
        change_type_counts[change_type] = change_type_counts.get(change_type, 0) + 1

    # Display in columns
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("##### 🔴 Severity Distribution")
        for severity in ['High', 'Moderate', 'Low']:
            count = severity_counts.get(severity, 0)
            if count > 0:
                if severity == 'High':
                    st.error(f"🔴 {severity}: {count}")
                elif severity == 'Moderate':
                    st.warning(f"🟡 {severity}: {count}")
                else:
                    st.success(f"🟢 {severity}: {count}")

    with col2:
        st.markdown("##### 📊 Change Types")
        sorted_changes = sorted(change_type_counts.items(), key=lambda x: x[1], reverse=True)
        for change_type, count in sorted_changes[:5]:  # Top 5 change types
            st.write(f"• {change_type}: {count}")

    # Drug binding impact summary
    high_impact = sum(1 for m in mutations if 'High' in m.get('druggability_impact', ''))
    if high_impact > 0:
        st.warning(f"{high_impact} mutations may significantly impact drug binding")

def display_mutation_discovery(wt_protein_sequence: str, chain_start: int = 1, max_mutations: int = 20) -> Optional[str]:
    """
    Mutation discovery interface for identifying known variants from public databases.
    Simplified design to match existing UI components.
    """
    if not wt_protein_sequence or len(wt_protein_sequence.strip()) < 10:
        st.warning("Please provide a valid protein sequence (minimum 10 amino acids)")
        return None

    # Discovery button
    discover_button = st.button(
        "Discover Mutations",
        key="discover_mutations_btn",
        type="primary",
        use_container_width=False
    )

    if discover_button:
        try:
            # Step 1: Protein identification
            with st.spinner("Identifying protein from public databases..."):
                protein_info = find_closest_protein_match(wt_protein_sequence)

            if "error" in protein_info:
                st.error(f"Protein identification failed: {protein_info['error']}")
                return None

  
            # Display protein match information
            similarity = protein_info.get('similarity', 0)
            user_seq_length = len(re.sub(r'[^A-Z]', '', wt_protein_sequence.upper()))
            full_protein_length = protein_info.get('length', 0)

            if similarity >= 0.95:
                st.success(f"Found: {protein_info['protein_name']} ({protein_info['accession']})")
            elif similarity >= 0.8:
                st.info(f"Found: {protein_info['protein_name']} ({protein_info['accession']}) - {similarity:.1%} similarity")
            else:
                st.warning(f"Found: {protein_info['protein_name']} ({protein_info['accession']}) - {similarity:.1%} similarity")

            # Show sequence mapping information
            if user_seq_length < full_protein_length:
                st.info(f"Your sequence ({user_seq_length} aa) maps to positions in the full protein ({full_protein_length} aa). Only mutations within your sequence region will be shown.")

            # Step 2: Variant retrieval
            with st.spinner("Retrieving known variants..."):
                mutations = discover_mutations_for_sequence(wt_protein_sequence, max_mutations, chain_start)

            if not mutations:
                st.warning("No mutations found for this protein.")
                st.info("This could indicate limited variant data in public databases or a unique protein sequence. You can still proceed with manually specified mutations.")
                return None

            # Success summary
            st.success(f"Found {len(mutations)} known mutations")

            # Results display
            df = format_mutations_for_display(mutations)

            # Color coding for severity
            def color_severity(val):
                if val == "High":
                    return "background-color: #ffebee"
                elif val == "Moderate":
                    return "background-color: #fff8e1"
                elif val == "Low":
                    return "background-color: #e8f5e8"
                return ""

            styled_df = df.style.applymap(color_severity, subset=['Severity'])
            st.dataframe(styled_df, use_container_width=True, hide_index=True)

            # Summary statistics
            severity_counts = {}
            for mut in mutations:
                severity = mut.get('severity', 'Unknown')
                severity_counts[severity] = severity_counts.get(severity, 0) + 1

            if severity_counts:
                st.markdown("**Severity Distribution:**")
                col1, col2, col3 = st.columns(3)
                with col1:
                    if severity_counts.get('High', 0) > 0:
                        st.metric("High", severity_counts['High'])
                with col2:
                    if severity_counts.get('Moderate', 0) > 0:
                        st.metric("Moderate", severity_counts['Moderate'])
                with col3:
                    if severity_counts.get('Low', 0) > 0:
                        st.metric("Low", severity_counts['Low'])

            # Quick copy section
            st.markdown("**Copy mutations for input:**")
            col1, col2 = st.columns(2)

            with col1:
                st.markdown("**All mutations:**")
                all_mutations = format_mutations_for_input(mutations[:10])
                st.code(all_mutations)

            with col2:
                st.markdown("**High-impact only:**")
                high_impact = [m for m in mutations if m.get('severity') in ['High', 'Moderate']]
                high_impact_formatted = format_mutations_for_input(high_impact[:10])
                st.code(high_impact_formatted)

            # Return formatted mutations for auto-population
            return format_mutations_for_input(mutations[:8])

        except Exception as e:
            st.error(f"Error during mutation discovery: {str(e)}")
            return None

    return None