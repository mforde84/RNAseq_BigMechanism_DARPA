# fast_count

usage

./fast_count_multi num_threads gtf_file bam_file(s) > output

Files included:

fast_count_multi - reports all counts and TPM, multithreading support

fast_count_deseq - reports gene counts in deseq compatible format, multithreading support

Benchmark: 8 core 1M annotations for 2Gb sorted reads ~30 seconds compared to ~28 minutes for bedtools multicov.

Requires bamtools API library at run time, and c++0x for compile.

git clone https://github.com/pezmaster31/bamtools

cd bamtools

mkdir build

cd build

cmake ..

make

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path to/lib

g++ -I bamtools/include/ -L bamtools/lib/ -o fast_count_multi fast_count_multi.cpp -lz -lbamtools -fpermissive -pthread -std=c++0x
