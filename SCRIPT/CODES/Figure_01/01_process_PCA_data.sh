#!/bin/bash
set -eu pipefail



# -------- ENCODE MWORDS data into 1mer(OHE) and Flexibility parameters ------

INPUTFILE="DATA/01_PCA/MWORDS_pca_input.fa"
OUTDIR="RESULTS/01_PCA_test/"
THREADS="2"


mkdir -p "$OUTDIR"
# ----------------------------------------------------------------------------

# -----------------------
# STEP 1. encode 1mer
# -----------------------
# python SCRIPTS/encode_kmer_zhou_update2.py -h
# usage: encode_kmer_zhou_update2.py [-h] input_fasta output_dir output_file k

# Encode DNA sequences from a FASTA file and save to an output file

# positional arguments:
#   input_fasta  Path to the input FASTA file
#   output_dir   Path to the output directory to save encoded sequences
#   output_file  Name of the output file to save encoded sequences
#   k            Length of k-mers to encode

# options:
#   -h, --help   show this help message and exit

python "SCRIPT/CODES/encode_kmer.py" "$INPUTFILE" "$OUTDIR" "MWORDS_pca_input" "1"


# ---------------------------
# STEP 2. encode flexibility
# ---------------------------
# DNAflexpy -h
# usage: DNAflexpy [-h] [--window-size WINDOW_SIZE] [--feature FEATURE] [--threads THREADS]
#                  [--outfile OUTFILE]
#                  input_file

# Process a multifasta file and calculate bendability of DNA sequence

# positional arguments:
#   input_file            Path to the input multifasta file

# options:
#   -h, --help            show this help message and exit
#   --window-size WINDOW_SIZE
#                         Size of the processing window [default: 10]
#   --feature FEATURE     Feature(s) to calculate ('DNaseI', 'NPP') [default: DNaseI]
#   --threads THREADS     Number of threads
#   --outfile OUTFILE     Output file name [optional]

for feat in DNaseI NPP twistDisp trx stiffness; 
    do 
        
        DNAflexpy --feature "$feat" --threads "$THREADS" --outfile "${OUTDIR}/MWORDS_pca_input_${feat}_0.tsv" "$INPUTFILE"
    done


echo "End of data preperation for PCA"