# tpmcalc

usage:

$ ./tpmcalc --help

$ ./tpmcalc {start-size-column} {annotation-column} {start-counts-column} {counts-file} {tpm-output} {logtpm-output}

$ ./tpmcalv 2 1 4 sample_data.txt sample_data_tpm.txt sample_data_log2tpm.txt

compile: g++ -o tpmcalc tpmcalc.cpp --std=c++1z

Calculates transcripts per million (TPM) and  log2 transformed TPM values from raw RNAseq counts. While it's possible to do these types of calculations in R, there are diminishing returns in terms of computational performance and simplicity as the number of samples increase. This program was written to support TPM calculations on very large data sets (e.g., ~12,000 TCGA transciptomes). The program is single threaded, and supports multiple sample count files. If processing a larger number of samples, it's recommended to run parallel instances of the program then trim and merge resulting files. 

