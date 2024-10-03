import sys
from Bio.Seq import Seq
from Bio import SeqIO

# Standard genetic code for translating DNA to protein
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',  # Isoleucine and Methionine
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',  # Threonine
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',  # Asparagine, Lysine
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  # Serine, Arginine
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',  # Leucine
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',  # Proline
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',  # Histidine, Glutamine
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',  # Arginine
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',  # Valine
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',  # Alanine
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',  # Aspartic acid, Glutamic acid
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',  # Glycine
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',  # Serine
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',  # Phenylalanine, Leucine
    'TAC':'Y', 'TAT':'Y',                        # Tyrosine
    'TGC':'C', 'TGT':'C',                        # Cysteine (C)
    'TGG':'W',                                   # Tryptophan
    'TAA':'_', 'TAG':'_', 'TGA':'_'              # Stop codons
}

# Function to translate a DNA sequence to a protein string
def translate_dna_to_protein(dna_sequence):
    protein = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if codon in CODON_TABLE:
            amino_acid = CODON_TABLE[codon]
            if amino_acid == '_':  # Stop codon found
                break
            protein.append(amino_acid)
    return ''.join(protein)

# Function to find all ORFs in a sequence for a given frame
def find_orfs_in_frame(sequence):
    orfs = []
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    for i in range(len(sequence) - 2):
        if sequence[i:i+3] == start_codon:
            for j in range(i, len(sequence) - 2, 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    orf_sequence = sequence[i:j+3]
                    protein = translate_dna_to_protein(orf_sequence)
                    if protein:  # Only add valid protein sequences
                        orfs.append(protein)
                    break
    return orfs

# Main function to find ORFs in all six reading frames
def find_all_orfs(dna_sequence):
    all_orfs = set()  # Use a set to avoid duplicate proteins
    
    # Find ORFs in the three forward frames
    for frame in range(3):
        orfs = find_orfs_in_frame(dna_sequence[frame:])
        all_orfs.update(orfs)
    
    # Get the reverse complement and find ORFs in the three reverse frames
    reverse_complement = str(Seq(dna_sequence).reverse_complement())
    for frame in range(3):
        orfs = find_orfs_in_frame(reverse_complement[frame:])
        all_orfs.update(orfs)
    
    return all_orfs

# Main function to read the FASTA file and output the proteins
def main(fasta_file, output_file, min_length):
    with open(output_file, 'a') as out_f:  # Open output file in append mode
        for record in SeqIO.parse(fasta_file, "fasta"):
            dna_sequence = str(record.seq)
            orfs = find_all_orfs(dna_sequence)
            for protein in orfs:
                if len(protein) >= min_length:
                    out_f.write(f">{record.id}\n{protein}\n")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python orf_finder.py <input_fasta> <min_length> <output_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]  # Input file in FASTA format
    min_length = int(sys.argv[2])  # Minimum length for ORF (in amino acids)
    output_file = sys.argv[3]  # Output file to store results
    main(fasta_file, output_file, min_length)
