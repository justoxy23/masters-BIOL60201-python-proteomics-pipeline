import argparse
import re

def parse_arguments():
    """
    Parse command-line arguments for the protein digestion script.

    Returns:
        argparse.Namespace: Parsed command-line arguments including:
            - input_file (str): Path to the input FASTA file.
            - enzyme (str): Enzyme code for digestion ('t', 'l', 'a', 'g').
            - missed_cleavages (int): Number of allowed missed cleavages (default: 0).
            - output_file (str): Path to the output FASTA file.
    """
    parser = argparse.ArgumentParser(
        description = "A Python script to digest protein sequences from a FASTA file using user-specified enzymes. The script supports missed cleavages and outputs the digested peptides to a file.",
        epilog = "Example usage: python3 group8_task2.py --input_file group8_task1_output.fasta --missed_cleavages 1 --enzyme t --output_file group8_task2_output.fasta")

    parser.add_argument('--input_file', required=True, help="Path to input fasta file."),
    parser.add_argument('--enzyme', required=True, choices=['t','l','a','g'], help="Enzyme used for digestion: t (trypsin), l (endoproteinase lys-c), a (endoproteinase arg-c), g (v8 proteinase)"),
    parser.add_argument('--missed_cleavages', type=validate_missed_cleavages, default = 0, help="Number of cleavages allowed (default:0, maximum:1)"),
    parser.add_argument('--output_file', required=True, help="Path to output fasta file.")

    return parser.parse_args()

def validate_missed_cleavages(value):
    """
    Validate the maximum allowed missed cleavages.

    Args:
        value (str): Value provided for missed cleavages.

    Returns:
        int: Validated integer value for missed cleavages.

    Raises:
        argparse.ArgumentTypeError: If the value is not within the range [0, 1].
    """
    int_value = int(value)
    if int_value < 0 or int_value > 1:
        raise argparse.ArgumentTypeError(
            f"Invalid value: {value}. Missed cleavages must be between 0 and 1 (inclusive).")
    return int_value

def fastaread(filename):
    """
    Reads a FASTA file and parses protein sequences.

    Args:
        filename (str): Path to the FASTA file.

    Returns:
        tuple: A tuple containing:
            - list: Order of protein identifiers in the file.
            - dict: Protein sequences keyed by their identifiers.

    Raises:
        ValueError: If the file is not found or contains no valid sequences.
    """   
    order = []
    sequences = {}
        
    try:
        with open(filename,"r") as fastafile:
            current_protein = None
            for line in fastafile:
                line = line.strip()
                if line.startswith(">"):
                    current_protein = line[1:].split()[0]
                    order.append(current_protein)
                    sequences[current_protein] = ""
                else:
                    sequences[current_protein] += line
    except FileNotFoundError:
        raise ValueError(f"Error: File '{filename}' not found.")
        
    if not order or not sequences:
        raise ValueError(f"Error: No valid protein sequences found in '{filename}'.")
        
    return order, sequences
        
def digest(sequence, enzyme, missed_cleavages=0):
    """
    Digest a protein sequence into peptides based on enzyme specificity and missed cleavages.

    Args:
        sequence (str): Protein sequence to digest.
        enzyme (str): Enzyme code specifying digestion rules ('t', 'l', 'a', 'g').
        missed_cleavages (int, optional): Number of missed cleavages allowed. Defaults to 0.

    Returns:
        list of tuple: A list of digested peptides, where each tuple contains:
            - str: Peptide sequence.
            - int: Number of missed cleavages in the peptide.

    Raises:
        ValueError: If an unsupported enzyme code is provided.
    """
    specificity = {
       "t" : r"([KR](?!P))", # Cuts at K or R unless the next residue is P
       "l" : r"([K](?!P))", # Cuts at Kunless the next residue is P
       "a" : r"([R](?!P))", # Cuts at R unless the next residue is P
       "g" : r"([E](?!P))", # Cuts at E unless the next residue is P
    }
    if enzyme not in specificity:
        raise ValueError(f"Unsupported enzyme: {enzyme}")
                   
    pattern = specificity[enzyme]
    sequence = re.sub(pattern, r"\1\n", sequence)
    peptides = sequence.split("\n")
    
    digested_peptides = []
    for i in range(len(peptides)):
        current_peptide = peptides[i]
        digested_peptides.append((current_peptide, 0))
        for j in range(1, missed_cleavages + 1):
            if i + j < len(peptides):
                current_peptide += peptides[i + j]
                digested_peptides.append((current_peptide, j))
        
    return digested_peptides
                         
def write_output(output_file, protein, peptides, enzyme):
    """
    Write digested peptides to a FASTA file.

    Args:
        output_file (str): Path to the output FASTA file.
        protein (str): Protein identifier.
        peptides (list of tuple): List of digested peptides with missed cleavages.
        enzyme (str): Enzyme code used for digestion.
    """
    enzyme_names = {
        "t": "Trypsin",
        "l": "Endoproteinase Lys-C",
        "a": "Endoproteinase Arg-C",
        "g": "V8 Proteinase (Glu-C)"
    }
    with open(output_file, "a") as file:           
        for idx, (peptide, missed) in enumerate(peptides, 1):
            file.write(f">{protein} peptide {idx} missed={missed} {enzyme_names[enzyme]}\n")
            file.write(f"{peptide}\n")
            
            
def main():
    """
    Main function to handle the script workflow:
    - Parse arguments.
    - Read input FASTA file.
    - Perform digestion.
    - Write output to file.
    """
    args = parse_arguments()
     
    try:
        order, seqs = fastaread(args.input_file)
    except ValueError as e:
        print(e)
        exit(1)
                         
    with open(args.output_file, "w") as file:
        pass
                         
    for protein in order:
        peptides = digest(seqs[protein], args.enzyme, args.missed_cleavages)
        write_output(args.output_file, protein, peptides, args.enzyme)
        
    print(f"Digestion completed. Results written to {args.output_file}")
        
if __name__ == "__main__":
    main()
                

# python task2-terenia.py --input_file task2-terenia-input.fasta --missed_cleavages 2 --enzyme t --output_file task2-terenia-output.fasta

