import sys

def flipComp(Text):
    # Compute the reverse complement of a DNA sequence.
    result = ""
    for letter in Text:
        if letter == 'A':
            result += 'T'
        elif letter == 'T':
            result += 'A'
        elif letter == 'G':
            result += 'C'
        else:
            result += 'G'
    flipped = result[::-1]  # Reverse the sequence
    return flipped

def read_fasta(filename):
    # Reads a single FASTA file and returns the sequence.
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

def find_genes(sequence, include_reverse=False):
    """Finds genes in reading frames (with or without reverse complements)."""
    forward_start_codon = 'ATG'
    forward_stop_codons = ['TAA', 'TAG', 'TGA']

    reverse_start_codon = 'CAT'  # Reverse complement of ATG
    reverse_stop_codons = ['TTA', 'CTA', 'TCA']  # Reverse complements of stop codons

    genes = []

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == forward_start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in forward_stop_codons:
                        gene = sequence[i:j+3]
                        genes.append(f"Forward: {gene}")
                        break

    # Check reverse reading frames if required
    if include_reverse:
        rev_sequence = flipComp(sequence)  # Get reverse complement using flipComp
        print(f"Reverse Complement Sequence: {rev_sequence[:50]}...")  # Print first 50 bases of reverse complement

        for frame in range(3):
            for i in range(frame, len(rev_sequence), 3):
                codon = rev_sequence[i:i+3]
                if codon == reverse_start_codon:  # Check for CAT in reverse complement
                    print(f"Found Reverse Start Codon (CAT) at position {i} in frame {frame}")  # Debug statement
                    for j in range(i+3, len(rev_sequence), 3):
                        stop_codon = rev_sequence[j:j+3]
                        if stop_codon in reverse_stop_codons:  # Look for TTA, CTA, TCA
                            gene = rev_sequence[i:j+3]
                            genes.append(f"Reverse: {gene}")
                            print(f"Found Reverse Gene from {i} to {j}")  # Debug statement
                            break

    return genes

def filter_genes_by_length(genes, min_length=100):
    """Filter genes by length (minimum codon length)."""
    filtered_genes = []
    for gene in genes:
        gene_sequence = gene.split(': ')[1]  # Extract the gene sequence
        if len(gene_sequence) >= min_length * 3:  # Apply length filtering
            filtered_genes.append(gene)
    return filtered_genes

def main():
    if len(sys.argv) != 2:
        print("Usage: python gene_finder_task2.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

    # Find genes in forward and reverse reading frames
    genes = find_genes(sequence, include_reverse=True)

    # Apply filtering by gene length
    filtered_genes = filter_genes_by_length(genes)

    output_file = 'output2_genes.txt'  # Updated to save to output2_genes.txt

    # Save both forward and reverse genes to the output file
    with open(output_file, 'w') as f:
        for gene in filtered_genes:
            if "Forward" in gene:
                f.write(f"Filtered Gene (Forward): {gene.split(': ')[1]}\n")
            elif "Reverse" in gene:
                f.write(f"Filtered Gene (Reverse): {gene.split(': ')[1]}\n")

    print(f"Total genes (forward and reverse) written: {len(filtered_genes)}")  # Debug statement

if __name__ == "__main__":
    main()

