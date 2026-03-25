def fastaread(filename):
    with open(filename) as fastafile:
        lines = fastafile.readlines()
    order = []
    seqs = {}
    for line in lines:
        if line.startswith('>'):
            items = line[1:].strip().split()
            name = items[0]
            number = items[2]
            order.append(items)
            dict_name = name + '_' + number
            seqs[dict_name] = []
        else:
            seqs[dict_name] = line.strip()
    return order, seqs


def pep2mass(seq, type):
    mono = {
        'A': 71.0371, 'C': 103.0092, 'D': 115.0269, 'E': 129.0426, 'F': 147.0684, 'G': 57.0215,
        'H': 137.0589, 'I': 113.0841, 'K': 128.0950, 'L': 113.0841, 'M': 131.0405, 'N': 114.0429,
        'P': 97.0528, 'Q': 128.0586, 'R': 156.1011, 'S': 87.0320, 'T': 101.0477, 'V': 99.0684,
        'W': 186.0793, 'Y': 163.0633, '*': 0.0
    }
    aver = {
        'A': 71.08, 'C': 103.14, 'D': 115.09, 'E': 129.12, 'F': 147.18, 'G': 57.05, 'H': 137.14, 'I': 113.16,
        'K': 128.17, 'L': 113.16, 'M': 131.19, 'N': 114.10, 'P': 97.12, 'Q': 128.13, 'R': 156.19, 'S': 87.08,
        'T': 101.10, 'V': 99.13, 'W': 186.21, 'Y': 163.18, '*': 0.0
    }
    waterM = 19.0106
    waterA = 19.0153
    non_standard = []
    reaction = False

    if type == 'mono':
        chosen = mono
        pepmass = waterM
    elif type == 'aver':
        chosen = aver
        pepmass = waterA

    for i in seq:
        if i not in chosen:
            reaction = True
            non_standard.append(i)
            continue
        else:
            pepmass += chosen[i]

    return pepmass, reaction, non_standard


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Convert peptide sequences to masses, read from a Fasta file')

    parser.add_argument('-i', '--input', help='Input fasta file', required=True)
    parser.add_argument('-o', '--output', help='Output fasta file', default="pepmasses.out")
    parser.add_argument('-t', '--type', help='Type of measuring, mono or aver', default='mono')
    parser.add_argument('-z', '--z_number', help='Peptide charge value', default=1, type=int)
    parser.add_argument(
        '-s', '--standard',
        help='Handling method of non-standard amino acids: (A) Ignore and report (B) Ignore without report '
             '(C) Do not output the mass-to-charge of this amino acid (D) Stop the program',
        default='A'
    )

    args = parser.parse_args()

    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Type of measuring: {args.type}")
    print(f"Peptide charge value: {args.z_number}")
    if args.standard == 'A':
        word = "(A) Ignore and report"
    elif args.standard == 'B':
        word = "(B) Ignore without reporting"
    elif args.standard == 'C':
        word = "(C) Do not output the mass-to-charge for the amino acid"
    elif args.standard == 'D':
        word = "(D) Stop the program"
    print(f"Handling method of non-standard amino acids: {word}")

    inputfile = args.input
    outputfile = args.output
    masstype = args.type
    ofile = open(outputfile, 'w')

    ion = 1.007
    z = args.z_number
    infor, seqs = fastaread(inputfile)
    masses = {}
    mass_to_change = {}

    for p_ID, seq in seqs.items():
        mass, reaction, non_standard = pep2mass(seq, masstype)

        if reaction:
            if args.standard == 'D':
                print("The sequence contains non-standard amino acids!!!\nThe program will EXIT automatically.")
                raise SystemExit
            elif args.standard == 'C':
                masses[p_ID] = ""
                mass_to_change[p_ID] = ""
                continue
            else:
                masses[p_ID] = mass
                if args.standard == 'A':
                    p_name = p_ID[:-2]
                    number = p_ID[-1]
                    print(f"Fragment No.{number} of {p_name} contains non-standard amino acids: {','.join(non_standard)}")
        else:
            masses[p_ID] = mass

        if z == 0:
            mass_to_change[p_ID] = mass
        else:
            mass_to_change[p_ID] = (mass + ion * z) / z

    ofile.write("Prot_name\tpeptide\tmass-to-charge\tz\tp\tsequence\n")
    for i in range(len(mass_to_change)):
        protein = infor[i][0]
        p_ID = list(seqs.keys())[i]
        w = protein + "\t"  # Prot_name
        w += infor[i][2] + "\t"  # Peptide
        if mass_to_change[p_ID] == "":
            w += mass_to_change[p_ID] + "\t"
        else:
            w += str(f"{mass_to_change[p_ID]:.3f}") + "\t"  # Mass-to-charge
        w += str(z) + "\t"  # Z
        w += infor[i][3][-1] + "\t"  # P
        w += list(seqs.values())[i] + "\n"  # Sequence
        ofile.write(w)

    ofile.close()
