from argparse import ArgumentParser
from Bio import SeqIO

def parse_fasta_file(file_path):
    """Read sequences from a FASTA file."""
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def find_orfs(sequence):
    """Find all open reading frames (ORFs) in a given sequence."""
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    # Check three different reading frames
    for frame in range(3):
        start_positions = []  # To store start codon positions
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i + 3]
            if codon == start_codon:
                start_positions.append(i)  # Found a start codon
            elif codon in stop_codons and start_positions:
                # Found a stop codon; create ORFs from all start positions
                while start_positions:
                    start_pos = start_positions.pop(0)
                    orf = sequence[start_pos:i + 3]
                    orfs.append((start_pos, i + 3, orf))

    return orfs

if __name__ == "__main__":
    parser = ArgumentParser(description="Find ORFs in a FASTA file.")
    parser.add_argument("-f", "--file", help="FASTA file", required=True)
    args = parser.parse_args()

    sequences = parse_fasta_file(args.file)
    
    for seq_id, sequence in sequences.items():
        print(seq_id)
        for start_pos, end_pos, orf in find_orfs(sequence):
            print(start_pos, end_pos, orf)
