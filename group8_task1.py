#!/usr/bin/env python
### Task 1: ORF finder 

import argparse                       ### Import for arguments 
import re, sys, os                    ### For additional stuff

### Define functions 


def fastaread(filename):             ### Reads the FASTA file and returns in a structured format: 
    
    """
    Reads a FASTA file and returns a tuple containing:
    - A list of sequence names (order)
    - A dictionary mapping sequence names to their sequences
    
    """
    fastafile = open( filename, "r")
    lines = fastafile.readlines()      
    fastafile.close()                 ### Read lines in one go into a list, 'lines'
      
    
    order = []                        ### List to keep the sequences in order 
    seqs = {}                         ### A dictionary thats stores all the sequences by name. 
    
    name = None 
    for line in lines:
        line = line.strip()               # Remove any whitespace and newlines 
        if line.startswith('>'):          # Header line 
            name = line [1:].split()[0]   # Take the header name 
            order.append(name)            # Add header to order list 
            seqs[name] = ''               # Start an empty sequence in the dictionary 
        elif name:
                seqs[name] += line.upper()               # Convert to uppercase 
        else:                                            # If no header is found, the file format is invalid
            print("This file is not in FASTA format")
            return None 
        
    print ('%d sequences read in' % len(order))  
    
    return (order, seqs)    # Return the list of names and sequence dictionary


def complement(s): 
    """
    Return the complementary sequence string.
    
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n' 
                     } 
    bases = list(s)                                                                          
    compseq = [basecomplement[base] for base in bases]    
    return ''.join(compseq)                            # Map each base to its complement  

def revcomplement(s):
    """Return the reverse complement sequence string."""
    s = s[::-1]                                         # Reverse the string and compute its complement
    s = complement(s)
    return s

def translate(dna):
    """
    Translate codons to amino acids for a list of ORFs.
    
    """
    gencode = { 
        "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", 
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T", 
        "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
        "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I", 
        "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", 
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", 
        "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A", 
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", 
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y", 
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", 
        "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
        "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"
    }

    amino_acids = []
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if 'N' in codon:
            print(f"Warning: Ignoring codon {codon} at position {i+1} due to ambiguous base N")
            return None                                          # Ignore ORFs with ambiguous bases
        if codon in gencode:
            amino_acids.append(gencode[codon])
        else:
            amino_acids.append('X')
            
    return ''.join(amino_acids)

#Define function to find ORFs in a DNA sequence
def find_orfs(dna, frame_offset, min_orf_length):
    """
    Find ORFs in a DNA sequence in all frames.
    Args: 
    - dna: The input DNA sequence (as a string)
    - frame_offset: The starting position for the reading frame (0, 1, or 2)
    - min_orf_length: Minimum length of an ORF 
    Returns: 
    - A list of tuples, each containing (start, end, sequence)
    
    """
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []                                           # List to store identified ORFs
    
    for i in range(frame_offset, len(dna), 3):          # Start scanning from the given offset
        if dna[i:i+3] == start_codon:
            for j in range(i, len(dna), 3):              
                codon = dna[j:j + 3]                       
                if codon in stop_codons: 
                    orf_sequence = dna[i:j+3]
                    amino_acid_sequence = translate(orf_sequence)
                    if amino_acid_sequence and len(amino_acid_sequence) >= min_orf_length:
                        orfs.append((i+1, j+3, amino_acid_sequence))
                    break
                    
    return orfs
    
def fasta_print(seq):         
    linelen = 80
    for pos in range(0,len(seq),linelen):
        print(seq[pos:pos+linelen])


# ## Test the function
# fastafile = "genome-small.fasta"
# result = fastaread(fastafile)

# ## Print the results 
# if result:
#     order, seqs = result
#     print('Order of sequences:', order)
#     print('Sequences dictionary:', seqs)

def argparse_min_length_ORF(value):
    """Validate that the input minimim ORF size is 50"""
    int_value = int(value)
    if int_value < 50:
        raise argparse.ArgumentTypeError(f"Minimum ORF size must be at least 50, but {int_value} is provided.")
    return int_value

def main():
    
    parser = argparse.ArgumentParser(description='Find ORFs in a nucleotide genome sequence')
    parser.add_argument("filename", help="Input FASTA file name.")
    parser.add_argument("-m", "--minOrfSize", help="the minimum length of an ORF in amino acids", type=argparse_min_length_ORF, default=50)
    parser.add_argument("-o", "--output", default="ORFs.fasta", help="Output file name")
    args = parser.parse_args()
                        
    inputfile = args.filename
    min_orf_length = args.minOrfSize
    outputfile = args.output
    
    order, seqs = fastaread(inputfile)
    
    ofile = open(outputfile, 'w')    
        
    for seq_name in order: 
        dna = seqs[seq_name]
        rev_dna = revcomplement(dna)
            
        orflist = []
        orfseqs = {}
            
        for frame in range(3):
            orfs = find_orfs(dna, frame, min_orf_length) 
            for idx, (start, end, sequence) in enumerate(orfs, start=1):
                orf_name = f"{seq_name}_F{frame + 1}_{idx:04d}"
                orflist.append((orf_name, frame+1, len(sequence), start))
                orfseqs[orf_name] = sequence
                
        for frame in range(3):
            orfs = find_orfs(rev_dna, frame, min_orf_length) 
            for idx, (start, end, sequence) in enumerate(orfs, start=1):
                orf_name = f"{seq_name}_F{frame + 4}_{idx:04d}"
                orflist.append((orf_name, frame+4, len(sequence), len(dna)-end+1))
                orfseqs[orf_name] = sequence

        # Write ORFs to the output file
        for orf_name, frame, length, start in orflist:
            header = f">{orf_name} frame:{frame} length:{length} start:{start}"
            ofile.write(f"{header}\n")
            ofile.write(f"{orfseqs[orf_name]}\n")

    ofile.close()
    print(f"ORFs have been written to {outputfile}")


if __name__ == "__main__":
    main()
else:
    print("run as module\n")