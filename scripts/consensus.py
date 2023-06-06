# This script is used to generate a consensus sequence from a FASTA file containing multiple DNA sequences.
# Usage python script.py input.fasta output.fasta
import sys


def generate_consensus(fasta_file):
    iupac_codes = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': 'AG',
        'Y': 'CT',
        'S': 'GC',
        'W': 'AT',
        'K': 'GT',
        'M': 'AC',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }

    # Read the FASTA file
    sequences = []
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                continue  # Skip the sequence identifier line
            sequences.append(line)

    # Check if all sequences have the same length
    sequence_length = len(sequences[0])
    if not all(len(seq) == sequence_length for seq in sequences):
        print("Error: Sequences are not aligned properly.")
        return None

    # Generate the consensus sequence
    consensus = ""
    for i in range(sequence_length):
        bases = [seq[i] for seq in sequences if seq[i] not in 'N-']
        base_counts = {base: bases.count(base) for base in set(bases)}
        if base_counts:
            if len(base_counts) == 1:
                consensus += next(iter(base_counts))
            else:
                ambiguity_code = ''.join(base for base in base_counts.keys() if base in iupac_codes)
                consensus += iupac_codes.get(ambiguity_code, 'N')
        else:
            consensus += 'N'  # If no valid base found, use 'N'

    return consensus


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python script.py input.fasta output.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    consensus_sequence = generate_consensus(input_file)
    if consensus_sequence:
        with open(output_file, 'w') as file:
            file.write(">Consensus\n")
            file.write(consensus_sequence)
        print("Consensus sequence generated and saved to", output_file)
