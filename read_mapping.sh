#!/bin/bash

# Note: look at the read_mapping2.sh for more flexibility with genome versions

### DESCRIPTION ###
# This script filters the raw fastq files and maps then maps them to the hg38 genome, creates bam files and indexes them. [fastq -> bam]
# The script takes 4 arguments:
#  i = input directory with the raw fastq files
#  r = run for which the analysis should be done
#  l = min read length to consider
#  o = output directory


while getopts i:r:l:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        r) run=${OPTARG};;
        l) len=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

mkdir $output/$run/
## filter out the reads that are shorter than 120 (amplicons1) or 105 (amplicons2) nt
#/usr/bin/parallel --jobs 6 --ungroup --link -- 'fname=$(basename {1} "R1_001.fastq.gz"); echo filtering $(basename $fname "_L001_"); /dir/bbmap/reformat.sh in={1} in2={2} out="/dir/fastq_filtered/'$run'/""$fname"R1_001.filtered.fastq.gz out2="/dir/'$run'/""$fname"R2_001.filtered.fastq.gz minlength='$len' 2> error_log.txt' ::: ${fws[@]} ::: ${rws[@]}


echo "finding paired reads from the filtered fastq files"
fws=$(find $input/$run"/" -name "*R1_001.filtered.fastq.gz" | sort)
rws=$(find $input/$run"/" -name "*R2_001.filtered.fastq.gz" | sort)

# map the reads to the genome, add read groups and filter out:
#   all reads where both paired reads didn't map (-f 0x2)
#   all reads which mapped to multiple locations (-F 0x104)
#   all reads whose mapping quality is <20 (-q 20)   
echo "mapping the reads to the hg38 with the BWA"
/usr/bin/parallel --jobs 6 --ungroup --link -- 'fname=$(basename {1} "_L001_R1_001.filtered.fastq.gz"); echo aligning $fname; id=$(zcat {1} | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g'); bwa mem -t 6 -v 3 "/hgXX.fa" {1} {2} -R $(echo "@RG\tID:$id\tSM:$fname\tLB:$id"_"$fname\tPL:ILLUMINA") | samtools view -b | samtools sort | samtools view -F 0x104 -f 0x2 -q 20 -b -o '$output'/'$run'/$fname.sorted.filtered.rg.bam' ::: ${fws[@]} ::: ${rws[@]} 

echo "indexing bam files"
bams=($(find $output/$run/"/" -name '*.sorted.filtered.rg.bam'))
parallel --jobs 6 --ungroup 'echo indexing {}; samtools index {} '$output'/'$run'/$(basename {}).bai' ::: ${bams[@]} 2> error_log.txt









