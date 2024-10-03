import sys

# Codon table mapping DNA codons to amino acids
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

# RBS sequence and upstream distance for Shine-Dalgarno sequence
RBS_SEQUENCE = "AGGAGG"
UPSTREAM_DISTANCE = 20

# Function to find reverse complement
def flipComp(sequence):
    """Compute the reverse complement of a DNA sequence, handling ambiguous nucleotide codes."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N'
    }
    return ''.join([complement.get(base, 'N') for base in sequence[::-1]])

# Function to read FASTA file
def read_fasta(filename):
    """Reads a single FASTA file and returns the sequence."""
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

# Function to translate DNA into a protein sequence
def translate_dna(dna_sequence):
    """Translates a DNA sequence into a protein sequence."""
    protein = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = CODON_TABLE.get(codon, "")
            if amino_acid == '_':  # Stop codon
                break
            protein += amino_acid
    return protein

# Function to check for RBS sequence
def find_rbs(sequence, start_index):
    """Check for RBS sequence within a specified upstream distance of the start codon."""
    upstream_region = sequence[max(0, start_index - UPSTREAM_DISTANCE):start_index]
    return RBS_SEQUENCE in upstream_region

# Function to find ORFs
def find_orfs(sequence, min_length=100):
    """Finds all ORFs in all six reading frames (both forward and reverse)."""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        if find_rbs(sequence, i) and (j - i + 3) // 3 >= min_length:
                            orf = sequence[i:j+3]
                            orfs.append(translate_dna(orf))
                        break

    # Check reverse reading frames
    rev_sequence = flipComp(sequence)
    for frame in range(3):
        for i in range(frame, len(rev_sequence), 3):
            codon = rev_sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(rev_sequence), 3):
                    stop_codon = rev_sequence[j:j+3]
                    if stop_codon in stop_codons:
                        if find_rbs(rev_sequence, i) and (j - i + 3) // 3 >= min_length:
                            orf = rev_sequence[i:j+3]
                            orfs.append(translate_dna(orf))
                        break

    return orfs

# Main function to execute the ORF finder
def main():
    if len(sys.argv) != 2:
        print("Usage: python orf_finder_with_rbs_filter.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

    # Find ORFs and filter them
    orfs = find_orfs(sequence)

    # Write the ORFs to output file
    output_file = f"./task6_filter_rbs_outputs/{fasta_file.split('/')[-1].replace('.fna', '_task6_filtered_orfs.txt')}"
    with open(output_file, 'w') as f:
        for orf in orfs:
            f.write(f"{orf}\n")

    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()



