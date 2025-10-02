# Biopython-Pipeline
GeneScope is a lightweight Streamlit-based bioinformatics tool that allows users to fetch DNA/RNA sequences from NCBI using accession IDs or FASTA uploads and perform core analyses in one place. It provides sequence statistics, GC/AT content, transcription/translation, codon usage, ORF detection, annotated CDS, BLASTn searches, and protein structure lookup with embedded 3D visualization. Designed with simplicity in mind, it uses only Biopython functions for computation and presents results with interactive visuals and tables. GeneScope is intended as an easy-to-use pipeline for students and researchers to replicate basic bioinformatics workflows on their local machines.

## Features
- Fetch sequence by NCBI accession or upload FASTA
- Sequence statistics (GC%, AT%, molecular weight, base counts, ambiguous %)
- Transcription & translation
- Optional ORF finding, codon usage, annotated CDS extraction
- Pairwise alignment with uploaded FASTA sequences
- BLASTn (top 5 hits)
- BLASTp vs PDB for structure lookup, embedded viewer
- Visualizations (base composition, ORF length distribution, codon usage)

## Installation
Clone this repo and install dependencies:

git clone https://github.com/yourusername/GeneScope.git
cd GeneScope
pip install -r requirements.txt

## Requirements
See `requirements.txt`

## License
MIT License
