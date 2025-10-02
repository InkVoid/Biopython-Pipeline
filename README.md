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

## Example Usage
<img width="1372" height="414" alt="Input   Settings" src="https://github.com/user-attachments/assets/45f97627-b27c-4595-bd7c-31c79cade79e" />

<img width="996" height="332" alt="Summary" src="https://github.com/user-attachments/assets/896619a1-5dc0-40de-948a-003b86a76576" />

<img width="980" height="661" alt="FASTA" src="https://github.com/user-attachments/assets/d89c11cd-5533-43d0-a8a1-ac513b4816f2" />

<img width="979" height="288" alt="Central Dogma" src="https://github.com/user-attachments/assets/411fe05f-3c41-432b-b0d3-642afe55dc14" />

<img width="1200" height="800" alt="base_composition" src="https://github.com/user-attachments/assets/f3da5e99-d17b-40c8-a955-a96a47297e01" />

<img width="1600" height="600" alt="codon_usage" src="https://github.com/user-attachments/assets/654190f6-9196-4b3c-b1fe-e6f5420f39e9" />

<img width="1200" height="800" alt="orf_histogram" src="https://github.com/user-attachments/assets/b36cc5e3-9b8a-4462-9c69-ec0cb79e9f48" />

## Installation
Clone this repo and install dependencies:

git clone https://github.com/yourusername/GeneScope.git
cd GeneScope
pip install -r requirements.txt

## Requirements
See `requirements.txt`

## License
MIT License
