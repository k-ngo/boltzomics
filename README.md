![BoltzOmics](https://github.com/k-ngo/boltzomics/blob/main/images/boltzomics.png?raw=true)
BoltzOmics is an interactive platform for deep learning-driven pharmacogenomic analysis. It leverages the Boltz-2 structure prediction model to assess how genetic mutations affect drug binding. Designed as a standalone Streamlit application, it offers an end-to-end solution for researchers in computational biology, pharmacology, and precision medicine, enabling the study of mutation effects on protein-ligand interactions.

![GUI](https://github.com/k-ngo/boltzomics/blob/main/images/gui.png?raw=true)
---

## Key Features

*   **State-of-the-Art Structure Prediction**: Utilizes Boltz-2 to generate high-quality 3D structures of protein-ligand complexes from sequence alone.
*   **Automated Mutation Discovery**: Queries the UniProt database to automatically find known mutations for your protein of interest. It integrates pre-computed SIFT and PolyPhen scores from the EBI Proteins API to provide insights into the functional impact of these variants.
*   **Manual Mutation Specification**: Allows users to manually input specific mutations of interest.
*   **Comprehensive Screening Modes**:
    *   *Protein-Drug Screening*: Test a panel of drugs against one or more proteins.
    *   *Mutation-Drug Screening*: Evaluate the effect of specific mutations on the binding of one or more drugs.
    *   *Structure-Only Mode*: Predict the structure of a protein without a ligand.
*   **Flexible and Customizable**: Supports multi-chain proteins, essential co-factors, post-translational modifications, and advanced modeling constraints.
*   **Interactive Visualization**: Includes an integrated 3D molecular viewer and a suite of plots to analyze binding affinity and structural quality metrics.
*   **User-Friendly Interface**: Implemented as a Streamlit application for easy local deployment and interaction.

---

## Setup and Installation

To set up and run BoltzOmics, follow these steps:

### Prerequisites

*   Python 3.12 or higher
*   Conda (Miniconda or Anaconda, see: https://www.anaconda.com/docs/getting-started/miniconda/install)
*   **NVIDIA GPU**: An NVIDIA GPU is highly recommended for optimal performance due to the computational intensity of the deep learning models.

### 1. Clone the Repository

First, clone this GitHub repository to your local machine by typing the following commands in the terminal:

```bash
git clone https://github.com/k-ngo/boltzomics.git
cd boltzomics
```

### 2. Create and Activate Conda Environment

In the terminal, create a new Conda environment and activate it:

```bash
conda create -n boltzomics python=3.12 -y
conda activate boltzomics
```

### 3. Install Dependencies

Install the required Python packages:

```bash
pip install -r requirements.txt
```

### 4. Run the Application

Once all dependencies are installed, you can launch the Streamlit application:

```bash
streamlit run boltzomics.py
```

This command will open the application in your web browser. If the browser does not automatically open or if running on a remote machine (e.g. through SSH connection), the Streamlit server will provide a URL (e.g., `http://localhost:8501` or an external IP address) that you can access from the web browser of your local machine to interact with the application.

---

## Usage

The application provides an intuitive interface for performing pharmacogenomic screenings.

1.  **Input Protein Sequence**: Provide your wild-type protein sequence or identifier.
2.  **Define Mutations**: Choose to either automatically discover known mutations from databases or manually specify mutations (e.g., `Y652A`, `F656A`).
3.  **Input Ligands**: Enter SMILES strings for the drugs you wish to screen.
4.  **Configure Parameters**: Adjust optional parameters for advanced modeling.
5.  **Run Prediction**: Initiate the structure prediction and screening process.
6.  **View Results**: Explore interactive plots, result tables, and 3D molecular visualizations of the predicted protein-ligand complexes.

---
