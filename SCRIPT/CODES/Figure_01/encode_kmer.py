from Bio import SeqIO
import itertools
import os
import argparse

def generate_kmers(k=1) -> list:
    """ Generate all possible k-mers of length k """
    return [''.join(p) for p in itertools.product('ACGT', repeat=k)]

def kmer_to_binary(kmer, kmers):
    """ Convert k-mer to a binary vector representation """
    vector = [0] * len(kmers)
    index = kmers.index(kmer)
    vector[index] = 1
    return vector
    
def sequence_to_feat(sequence, k):
    """ Convert a sequence to a binary vector based on k-mers """
    L = len(sequence)
    kmers = generate_kmers(k)
    binary_vector = []

    for i in range(L - k + 1):
        kmer = sequence[i:i + k]
        kmb = kmer_to_binary(kmer, kmers)
        binary_vector.extend(kmb)
        
    return binary_vector

def encode_sequences(sequences, k):
    """ Encode all sequences from a dictionary """
    encoded_sequences = {}
    for seqid, sequence in sequences.items():
        encoded_sequences[seqid] = sequence_to_feat(sequence, k)
    return encoded_sequences

def save_features(encoded_sequences, outdir, outfilepath, k):
    """ Save encoded sequences to an output file """
    outfile = os.path.join(outdir, f"{outfilepath}_{k}mer.tsv")
    
    with open(outfile, 'w') as f:
        for seqid, encoded_seq in encoded_sequences.items():
            encoded_str = '\t'.join(map(str, encoded_seq))
            f.write(f"{seqid}\t{encoded_str}\n")    

def read_fasta(file_path):
    """ Read sequences and IDs from a FASTA file using Biopython """
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Encode DNA sequences from a FASTA file and save to an output file")
    parser.add_argument("input_fasta", help="Path to the input FASTA file")
    parser.add_argument("output_dir", help="Path to the output directory to save encoded sequences")
    parser.add_argument("output_file", help="Name of the output file to save encoded sequences")
    parser.add_argument("k", type=int, help="Length of k-mers to encode")
    
    args = parser.parse_args()
    
    # Read sequences from the input FASTA file
    sequences = read_fasta(args.input_fasta)
    
    # Encode the sequences
    encoded_sequences = encode_sequences(sequences, args.k)
    
    # Save the encoded sequences to the output file
    save_features(encoded_sequences, args.output_dir, args.output_file, args.k)
    
if __name__ == "__main__":
    main()