#!/bin/bash

while getopts i: flag
do
    case "${flag}" in
        i) in_dir=${OPTARG};;
    esac
done


# setting the parameters
PICARD="java -jar /dir/picard.jar"


for dir in $(ls -1 $in_dir | grep -e "pattern"); do

  echo "renaming the samples in the directory:" $dir :

  for vcf in $(ls -d $in_dir/$dir/*); do
  
     echo renaming $(basename $vcf) 
     sample=$(basename $vcf .vcf) 
     $PICARD RenameSampleInVcf -I $vcf -O $in_dir/$dir/$(basename $vcf .vcf).renamed.vcf --NEW_SAMPLE_NAME $sample 2> $in_dir/renaming_error_log.txt
     mv $in_dir/$dir/$(basename $vcf .vcf).renamed.vcf $vcf
     bgzip $vcf
     tabix $vcf.gz
    
  done
  
  
  echo "merging the samples in the directory:" $dir : 
  
  vcfs=$(ls -d $in_dir/$dir/*.vcf.gz)
  bcftools merge $vcfs -o $in_dir/$dir.vcf 2> $in_dir/merging_error_log.txt
  
done






