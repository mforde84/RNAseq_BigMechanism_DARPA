#!/bin/bash
#usage ./run_pipeline.sh <genome director> <fastq directory>
#dependancies are based upon Ubuntu LTS 14.04. Different environments 
#will likely require additional modifications.

######
#env settings
######
sudo sh -c "echo 'kernel.shmmax = 31000000000' >> /etc/sysctl.conf";
sudo sh -c "echo 'kernel.shmall = 31000000000' >> /etc/sysctl.conf";
sudo /sbin/sysctl -p;
export threads=$(grep -c ^processor /proc/cpuinfo);
export threadstwo=$(($threads*2));
export org_dir="$(pwd)";
export genome_dir=$1;
export fastq_dir=$2;

######
#repo dependancies
######
sudo apt-get update;
sudo apt-get install build-essential cmake samtools -y;

######
#source dependancies
######
cd bamtools;
mkdir build;
cd build;
cmake ..;
make;
cd ../../mmr;
./configure;
make;
cd ../tpmcalc;
g++ -o tpmcalc tpmcalc.cpp --std=c++1z;
chmod +x tpmcalc;
cd ../fast_count;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$org_dir"/bamtools/lib;
g++ -I "$org_dir"/bamtools/include/ -L "$org_dir"/bamtools/lib/ -o fast_count_multi fast_count_multi.cpp -lz -lbamtools -fpermissive -pthread -std=c++0x;
chmod +x fast_count_multi;
export PATH=$PATH:"$org_dir"/STAR/bin/Linux_x86_64/:"$org_dir"/mmr/:"$org_dir"/fast_count/:"$org_dir"/tpmcalc/;

######
#add swap (need atleast 32 gb of memory to index hg20)
######
cd "$org_dir";
sudo dd if=/dev/zero of=swapper bs=1G count=32;
sudo mkswap swapper;
sudo swapon swapper;

######
#prepare assembly
######
cd "$genome_dir";
gunzip *.gz;
mkdir 2pass;
ln *.fa 2pass/genome.fa
ln *.gtf 2pass/genes.gtf

######
#generate genome index for first pass for reads length 50bp
######
STAR --runThreadN "$threadstwo" --runMode genomeGenerate --genomeDir . --genomeFastaFiles *.fa --sjdbGTFfile *.gtf --sjdbOverhang 49;

######
#generate alignments for first pass
######
cd "$fastq_dir";
gunzip *.gz;
for z in *.fastq; do
	STAR --runMode alignReads --outFileNamePrefix "$z".1pass. --runThreadN "$threadstwo" --genomeDir "$genome_dir" --genomeLoad LoadAndKeep --readFilesIn "$z" --outSAMtype BAM Unsorted --outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000  --alignMatesGapMax 1000000;
done;
STAR --genomeLoad Remove --outFileNamePrefix "$genome_dir"/genome.remove. --genomeDir "$genome_dir";

######
#combine novel splice junctions
######
for f in *.SJ.out.tab; do
 cat "$f" >> combo.SJ.out.tab;
done;

######
#generate genome index for second pass for reads length 50bp
######
cd "$genome_dir"/2pass;
STAR --runThreadN "$threadstwo" --runMode genomeGenerate --genomeDir "$genome_dir"/2pass --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbFileChrStartEnd combo.SJ.out.tab --sjdbOverhang 49;

######
#generate alignments for 2nd pass
######
cd "$fastq_dir";
for z in *.fastq; do
 STAR --runMode alignReads --outFileNamePrefix "$z".2pass. --runThreadN "$threadstwo" --genomeDir "$genome_dir"/2pass --genomeLoad LoadAndKeep --readFilesIn "$z" --outSAMtype BAM Unsorted --outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000  --alignMatesGapMax 1000000;
done;
STAR --genomeLoad Remove --outFileNamePrefix "$genome_dir"/2pass/genome.remove. --genomeDir "$genome_dir"/2pass;

######
#sort 2nd pass alignments
######
find fastq/ -name "*2pass*Aligned.out.bam" | xargs -n 1 -P "$threadstwo" -iFILES sh -c 'samtools sort -n FILES FILES.sort;';

######
#resolve multimapped reads
######
for f in *.sort.bam; do
 mmr -o "$f".mmr.bam -t "$threadstwo" -S -b "$f";
done;

######
#sort mmr alignments
######
find . -name "*.mmr.bam" | xargs -n 1 -P "$threadstwo" -iFILES sh -c 'samtools sort FILES FILES.sort;';

######
#count features and intergenic, fast_count is a bit heavy on memory, and 
#will use anywhere from 3-10gb of memory per process so don't over do it
######
find . -name "*mmr.bam.sort.bam" | xargs -n 1 -P 4 -iFILES sh -c 'fast_count_multi 1 $genome_dir/2pass/genes.gtf FILES > FILES.count.feature.txt; fast_count_multi 1 $org_dir/ident-intergenic.gtf FILES > FILES.count.intergenic.txt;';

######
#filter by gene feature
######
for f in *count*feature*; do
	python gene_puller.py "$f".pulls.txt;
	cut -f 10 $f.pulls.txt > $f.hits.txt;
done;

######
#extract feature annotation information
######
extract="$(ls *.pulls.txt | head -n 1)";
cut -f 4-5 "$extract" > feat.gen.loc;
cut -f 9 "$extract" | sed 's/^gene_id.\"//g' | sed 's/\".*//g' > feat.anno.loc
paste feat.anno.loc feat.gen.log | pr -t > feat.anno;

######
#extract intergenic annotation information
######
extract="$(ls *.intergenic.txt | head -n 1)";
cut -f 4-5 "$extract" > inter.gen.anno;
cut -f 9 "$extract" > inter.anno.loc;
paste inter.anno.loc inter.gen.log | pr -t > inter.anno;
for f in *intergenic.txt; do
	cut -f 10 $f.pulls.txt > $f.hits.txt;
done;

######
#merge raw counts
######
mkdir "$org_dir"/output;
paste feat.anno $(ls *feature*hits*) | pr -t > "$org_dir"/output/raw.feature.counts.txt;
paste inter.anno $(ls *intergenic*hits*) | pr -t > "$org_dir"/output/raw.intergenic.counts.txt;

######
#calc tpm
######
cd "$org_dir"/output/;
find . -name "*counts.txt" | xargs -n 1 -P 2 -iFILES sh -c 'tpmcalc 2 1 4 FILES FILES.tpm.txt FILES.log2tpm.txt;';

exit 0;
