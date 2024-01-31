#!/bin/bash

### DESCRIPTION ###
# This script maps the paired reads (fastq) to the reference genome and optionally removes the duplicated reads. (Indexing included)
# The script takes 3 arguments:
#  v = reference genome version
#  i = directory with the input (fastq) files
#  d = whether to deduplicate reads or not (by default it is set to true; any other value will turn it off)


# set the default value for the read deduplication
deduplicate=${d:-false}

while getopts v:i:d: flag
do
    case "${flag}" in
        v) version=${OPTARG};;
        i) in_dir=${OPTARG};;
        d) deduplicate=${OPTARG};;
    esac
done

echo "Genome version: $version"
reference="/dir/$version/$version.fa"



for dir in $(ls -1 $in_dir | grep -e "serija"); do
  echo processing $dir :
  
  out_dir="/dir/$version/$dir/"
  
  echo "finding reads"
  fws=$(find /$dir/ -name "*R1_001.fastq.gz" | sort)
  rws=$(find /$dir/ -name "*R2_001.fastq.gz" | sort)
  
  echo "mapping the reads with the BWA"
  parallel --jobs 6 --ungroup --link -- 'header=$(zcat {1} | head -n 1); id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g'); sm=$(basename {1} _L001_R1_001.fastq.gz); echo aligning $sm; bwa mem -t 6 -v 1 -R $(echo "@RG\tID:$id"_"$sm\tSM:$sm\tLB:$dir\tPL:ILLUMINA") '$reference' {1} {2} | samtools view -b | samtools sort | samtools view -q 20 -b -o '$out_dir'/$sm.sorted.filtered.bam' ::: ${fws[@]} ::: ${rws[@]} 2> mapping_error_log.txt
  
  
  if [ $deduplicate == true ]
  then
  
    ### marking duplicated reads ###
    echo "marking duplicated reads in the bam files..."
    bams=($(find $out_dir -name '*.sorted.filtered.bam'))
    parallel --jobs 6 --ungroup 'echo processing $(basename {}); java -jar /picard.jar MarkDuplicates -I {} -M '$out_dir'$(basename {} .sorted.filtered.bam)_metrics.txt -O '$out_dir'$(basename {} .bam).marked_dup.bam' ::: ${bams[@]} 2> marking_dups_error_log.txt
  
    echo "indexing bam files..."
    bams_dedup=($(find $out_dir -name '*.sorted.filtered.marked_dup.bam'))
    parallel --jobs 6 --ungroup 'echo processing $(basename {}); samtools index {}' ::: ${bams_dedup[@]} 2> indexing_error_log.txt
  
  else

    echo "indexing bam files..."
    bams=($(find $out_dir -name '*.sorted.filtered.bam'))
    parallel --jobs 6 --ungroup 'echo indexing {}; samtools index {}' ::: ${bams[@]} 2> indexing_error_log.txt
  
  fi

  
done






    