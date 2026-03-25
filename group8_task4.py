#!/usr/bin/env python3

import argparse
import sys
import csv
import numpy as np

def load_peptide_data(input_file):
    """
    Load peptide mass data from input file.
    
    Args:
        input_file (str): Path to input file with peptide mass data
    
    Returns:
        list: List of dictionaries containing peptide information
    """
    peptides = []
    try:
        with open(input_file, 'r') as f:
            # Skip header line
            next(f)
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                peptide_info = {
                    'prot_name': row[0],
                    'peptide': int(row[1]),
                    'mass_to_charge': float(row[2]),
                    'z': int(row[3]),
                    'missed_cleavages': int(row[4]),
                    'sequence': row[5]
                }
                peptides.append(peptide_info)
    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    
    return peptides

def mode1_peptide_count(peptides, min_mz, max_mz):
    """
    Count peptides within a specified m/z range.
    
    Args:
        peptides (list): List of peptide dictionaries
        min_mz (float): Minimum m/z value
        max_mz (float): Maximum m/z value
    
    Returns:
        int: Number of peptides in the specified range
    """
#     return sum(1 for p in peptides if min_mz <= p['mass_to_charge'] <= max_mz)

    protein_peptide_counts = {}

    # Filter peptides within the m/z range
    for peptide in peptides:
        if min_mz <= peptide['mass_to_charge'] <= max_mz:
            protein = peptide['prot_name']
            if protein not in protein_peptide_counts:
                protein_peptide_counts[protein] = 0
            protein_peptide_counts[protein] += 1

    return protein_peptide_counts

def mode2_peptide_histogram(peptides, bin_size):
    """
    Create a histogram of peptides based on user-defined bin size.
    
    Args:
        peptides (list): List of peptide dictionaries
        bin_size (float): Size of histogram bins
    
    Returns:
        dict: Histogram of peptide counts per bin
    """
    masses = [p['mass_to_charge'] for p in peptides]
    min_mass = min(masses)
    max_mass = max(masses)
    
    bins = np.arange(min_mass, max_mass + bin_size, bin_size)
    hist, _ = np.histogram(masses, bins=bins)
    
    return {f"{bins[i]:.2f}-{bins[i+1]:.2f}": count for i, count in enumerate(hist)}

def mode3_sliding_window(peptides, window_size, step_size):
    """
    Analyze peptides using a sliding window approach.
    
    Args:
        peptides (list): List of peptide dictionaries
        window_size (float): Size of the sliding window
        step_size (float): Step size between windows
    
    Returns:
        list: Peptide counts for each window
    """
    masses = [p['mass_to_charge'] for p in peptides]
    min_mass = min(masses)
    max_mass = max(masses)
    
    window_stats = []
    current_start = min_mass
    
    while current_start <= max_mass:
        window_end = current_start + window_size
        count = sum(1 for m in masses if current_start <= m < window_end)
        window_stats.append({
            'start': current_start,
            'end': window_end,
            'peptide_count': count
        })
        current_start += step_size
    
    return window_stats

def mode4_unique_protein_peptides(peptides, mass_accuracy):
    """
    Identify peptides unique to specific proteins.
    
    Args:
        peptides (list): List of peptide dictionaries
        mass_accuracy (float): Mass accuracy threshold
    
    Returns:
        dict: Proteins with unique peptides
    """
    # Group peptides by protein
    protein_peptides = {}
    for p in peptides:
        if p['prot_name'] not in protein_peptides:
            protein_peptides[p['prot_name']] = []
        protein_peptides[p['prot_name']].append(p['mass_to_charge'])
    
    # Check for unique peptides
    unique_proteins = {}
    for prot, masses in protein_peptides.items():
        unique_peptides = []
        for mass in masses:
            # Check if this mass is unique across all proteins within mass accuracy
            is_unique = all(
                not any(abs(mass - other_mass) <= mass_accuracy 
                        for other_mass in other_masses)
                for other_prot, other_masses in protein_peptides.items() 
                if other_prot != prot
            )
            if is_unique:
                unique_peptides.append(mass)
        
        if unique_peptides:
            unique_proteins[prot] = unique_peptides
    
    return unique_proteins

def main():
    parser = argparse.ArgumentParser(description='Peptide Mass Statistics Calculator')
    parser.add_argument('-i', '--input', required=True, 
                        help='Input file with peptide mass data')
    parser.add_argument('-o', '--output', 
                        help='Output file for results')
    parser.add_argument('-m', '--mode', type=int, choices=[1, 2, 3, 4], required=True,
                        help='Analysis mode: 1=range count, 2=histogram, 3=sliding window, 4=unique peptides')
    # Mode-specific arguments
    parser.add_argument('--min-mz', type=float, default=1000, 
                        help='Minimum m/z for mode 1 (default: 1000)')
    parser.add_argument('--max-mz', type=float, default=1500, 
                        help='Maximum m/z for mode 1 (default: 1500)')
    parser.add_argument('--bin-size', type=float, default=100, 
                        help='Bin size for mode 2 histogram (default: 100)')
    parser.add_argument('--window-size', type=float, default=200, 
                        help='Window size for mode 3 (default: 200)')
    parser.add_argument('--step-size', type=float, default=50, 
                        help='Step size for mode 3 sliding window (default: 50)')
    parser.add_argument('--mass-accuracy', type=float, default=0.1, 
                        help='Mass accuracy for mode 4 unique peptide detection (default: 0.1)')
    
    args = parser.parse_args()
    
    # Load peptide data
    peptides = load_peptide_data(args.input)
    
#     # Perform analysis based on selected mode
#     if args.mode == 1:
#         result = mode1_peptide_count(peptides, args.min_mz, args.max_mz)
#         output_text = f"Peptides in m/z range {args.min_mz}-{args.max_mz}: {result}"

    if args.mode == 1:
        # Get proteins and their respective peptide counts within the specified m/z range
        protein_peptide_counts = mode1_peptide_count(peptides, args.min_mz, args.max_mz)
        # Format the output text
        output_text = f"Proteins with peptides in m/z range {args.min_mz}-{args.max_mz}:\n"
        for protein, count in protein_peptide_counts.items():
            output_text += f"{protein}: {count} peptides\n"
    
    elif args.mode == 2:
        result = mode2_peptide_histogram(peptides, args.bin_size)
        output_text = "Peptide Histogram (bin: mass range, count):\n" + \
                      "\n".join(f"{bin}: {count}" for bin, count in result.items())
    
    elif args.mode == 3:
        result = mode3_sliding_window(peptides, args.window_size, args.step_size)
        output_text = "Sliding Window Analysis:\n" + \
                      "\n".join(f"Window {r['start']:.2f}-{r['end']:.2f}: {r['peptide_count']} peptides" 
                                for r in result)
    
    elif args.mode == 4:
        result = mode4_unique_protein_peptides(peptides, args.mass_accuracy)
        output_text = "Proteins with Unique Peptides:\n" + \
                      "\n".join(f"{prot}: {len(masses)} unique peptides" 
                                for prot, masses in result.items())
    
    # Output results
    if args.output:
        try:
            with open(args.output, 'w') as f:
                f.write(output_text)
            print(f"Results written to {args.output}")
        except Exception as e:
            print(f"Error writing to output file: {e}")
    else:
        print(output_text)

if __name__ == "__main__":
    main()