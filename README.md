# BIOL60201 Proteomics Pipeline (Python)

This repository contains a Python-based workflow developed to perform an end-to-end **proteomics analysis**, progressing from genomic sequences to peptide and ion-level outputs.

The pipeline is structured across four tasks, each representing a key stage in a typical proteomics workflow.

---

## Project overview

This project implements a simplified computational proteomics pipeline, including:

- Open Reading Frame (ORF) identification and translation
- Enzymatic digestion of protein sequences
- Peptide mass calculation
- Ion-level statistics and analysis

The workflow simulates key steps used in mass spectrometry-based proteomics and demonstrates how raw sequence data can be transformed into analytically meaningful outputs.

---

## Pipeline structure

The workflow is divided into four main components:

### Task 1: ORF identification and translation
- Identifies open reading frames from input nucleotide sequences  
- Translates nucleotide sequences into protein sequences  
- Outputs predicted protein sequences for downstream analysis  

### Task 2: Enzymatic digestion
- Simulates protease digestion (e.g. trypsin)  
- Splits protein sequences into peptides based on cleavage rules  
- Generates peptide lists for each protein  

### Task 3: Peptide mass calculation
- Computes molecular weights of peptides  
- May include m/z calculations depending on implementation  
- Prepares peptides for mass spectrometry interpretation  

### Task 4: Ion statistics and analysis
- Analyses peptide/ion distributions  
- Summarises properties relevant to detection in MS experiments  
- Produces outputs that support interpretation of peptide profiles  

---

## Input

The pipeline expects sequence data such as:

- DNA sequences (for ORF identification)
- Protein sequences (if starting from translated data)

Input formats may include standard text or FASTA-like structures depending on the script implementation.

---

## Output

The workflow generates outputs such as:

- Predicted protein sequences  
- Digested peptide lists  
- Peptide mass tables  
- Ion-level summaries and statistics  

These outputs collectively represent the transformation from genomic input to proteomics-ready data.

---

## Authors

Task 1
Johanna Olaniyi  
MSc Bioinformatics and Systems Biology  
University of Manchester

Task 2
Terenia Lee Leh Yi  
MSc Bioinformatics and Systems Biology  
University of Manchester

Task 3
Jia Yi Lyu  
MSc Bioinformatics and Systems Biology  
University of Manchester

Task 4
Sam Wheeler  
MSc Bioinformatics and Systems Biology  
University of Manchester

