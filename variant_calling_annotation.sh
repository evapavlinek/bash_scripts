#!/bin/bash

### DESCRIPTION ###
# This script calls the variants (either all or in targeted regions) and annotates them with the data from dbSnp, Clinvar, phastCons and dbNSFP databases. It creates a temporary directory in the output directory which will be used for 
# The script takes 5 arguments:
#  v = reference genome version
#  i = directory with the input (bam) files
#  r = path to the bed file with the coordinates of the targeted regions (by default it is set to false)
#  k = whether keep all files (temporary files created after certain steps) (by default it is set to false; any other value will turn it on)
#  o = output directory

# set the default value for deleting the temporary files
regions=${r:-false}
keep=${k:-false}

while getopts v:i:r:k:o: flag
do
    case "${flag}" in
        v) version=${OPTARG};;
        i) in_dir=${OPTARG};;
        r) regions=${OPTARG};;
        k) keep=${OPTARG};;
        o) out_dir=${OPTARG};;
    esac
done

# create the temporary directory
mkdir $out_dir/tmp
temp=$out_dir/tmp/


echo "Genome version: $version"
echo "Bed file with the regions of interest: $regions"

if [ $version == "hg19" ] 
then
  
  ### variant calling
  
  bams=($(find $in_dir -name '*.sorted.filtered.bam'))
  
  if [ $regions == "false" ] 
  then
    
    # this command calls all variants 
    echo "Calling variants..."
    parallel --jobs 6 --ungroup 'echo processing $(basename {}); freebayes -f /hg19.fa {} > '$temp'$(basename {} .bam).vcf' ::: ${bams[@]} 2> error_log.txt
    
  else
    
    # this command calls variants only in regions of interest/gene panel (provided BED file)
    echo "Calling variants only in the targeted regions..."
    parallel --jobs 6 --ungroup 'echo processing $(basename {}); freebayes -f /hg19.fa -t '$regions' {} > '$temp'$(basename {} .bam).targeted_regions.vcf' ::: ${bams[@]} 2> error_log.txt
  
  fi
  
  ### variant annotation
  
  # SnpEff: add ANN field 
  echo "annotating variants with the SnpEff..."
  vcfs1=($(find $temp -name '*.sorted.filtered.targeted_regions.vcf'))
  parallel --jobs 6 --ungroup 'echo processing $(basename {}); java -Xmx8g -jar /snpEff.jar -v -noStats hg19 {} > '$temp'$(basename {} .vcf).ann.vcf' ::: ${vcfs1[@]} 2> error_log.txt

  
  #SnpSift: add information from various databases
  echo "SnpSift: annotating variants with the dbSnp, Clinvar and phastCons databases..."
  vcfs2=($(find $temp -name '*.sorted.filtered.targeted_regions.ann.vcf'))
  parallel --jobs 6 --ungroup 'echo processing $(basename {}); bcftools sort {} | java -jar /SnpSift.jar annotate /hg19/dbSnp.vcf.gz | java -jar SnpSift.jar annotate /clinvar.vcf.gz | java -jar /SnpSift.jar phastCons /phastCons/ - > '$temp'$(basename {} .vcf).dbSnp.Clinvar.phastCons.vcf' ::: ${vcfs2[@]} 2> error_log.txt
  
  
  # this part is separated because the dbNSFP annotation won't work in the pipeline 
  echo "SnpSift: annotating variants with the dbNSFP scores and adding variant type and zygosity information..."
  vcfs3=($(find $temp -name '*.dbSnp.Clinvar.phastCons.vcf'))
  parallel --jobs 6 --ungroup 'echo processing $(basename {}); java -jar /SnpSift.jar dbnsfp -db /hg19/dbNSFP.txt.gz {} | java -jar /SnpSift.jar varType - > '$'$(basename {} .vcf).dbNSFP.vcf' ::: ${vcfs3[@]} 2> error_log.txt
  
  
  # deleting temporary files
  if [ $keep == "false" ] 
  then
    
    echo "Deleting temporary files..."
    rm -rf $temp   
    
  fi


elif [ $version == "hg38" ]
then
  
  ### variant calling
  
  bams=($(find $in_dir -name '*.sorted.filtered.bam'))
  
  if [ $regions == "false" ] 
  then
    
    # this command calls all variants 
    echo "Calling variants..."
    parallel --jobs 6 --ungroup 'echo processing $(basename {}); freebayes -f /hg38.fa {} > '$temp'$(basename {} .bam).vcf' ::: ${bams[@]} 2> error_log.txt
    
  else
    
    # this command calls variants only in regions of interest/gene panel (provided BED file)
    echo "Calling variants only in the targeted regions..."
    parallel --jobs 6 --ungroup 'echo processing $(basename {}); freebayes -f /hg38.fa -t '$regions' {} > '$temp'$(basename {} .bam).targeted_regions.vcf' ::: ${bams[@]} 2> error_log.txt
  
  fi
  
  ### variant annotation
  
  # SnpEff: add ANN field 
  echo "annotating variants with the SnpEff..."
  vcfs1=($(find $temp -name '*.sorted.filtered.targeted_regions.vcf'))
  parallel --jobs 6 --ungroup 'echo processing $(basename {}); java -Xmx8g -jar /snpEff.jar -v -noStats hg38 {} > '$temp'$(basename {} .vcf).ann.vcf' ::: ${vcfs1[@]} 2> error_log.txt

  
  #SnpSift: add information from various databases
  echo "SnpSift: annotating variants with the dbSnp, Clinvar and phastCons databases..."
  vcfs2=($(find $temp -name '*.sorted.filtered.targeted_regions.ann.vcf'))
  parallel --jobs 6 --ungroup 'echo processing $(basename {}); bcftools sort {} | java -jar /SnpSift.jar annotate /dbSNP/00-All.vcf.gz | java -jar /SnpSift.jar annotate /ClinVar/clinvar.vcf.gz | java -jar /SnpSift.jar phastCons /hg38/phastCons/ - > '$temp'$(basename {} .vcf).dbSnp.Clinvar.phastCons.vcf' ::: ${vcfs2[@]} 2> error_log.txt
  
  
  # this part is separated because the dbNSFP annotation won't work in the pipeline 
  echo "SnpSift: annotating variants with the dbNSFP scores and adding variant type and zygosity information..."
  vcfs3=($(find $temp -name '*.dbSnp.Clinvar.phastCons.vcf'))
  parallel --jobs 6 --ungroup 'echo processing $(basename {}); java -jar SnpSift.jar dbnsfp -db /hg38/dbNSFP/dbNSFP.txt.gz {} | java -jar /SnpSift.jar varType - > '$out_dir'$(basename {} .vcf).dbNSFP.vcf' ::: ${vcfs3[@]} 2> error_log.txt
  
  
  # deleting temporary files
  if [ $keep == "false" ] 
  then
    
    echo "Deleting temporary files..."
    rm -rf $temp   
    
  fi

fi
