#!/bin/bash

### DESCRIPTION ###
# This script performs the joint variant calling (either on all variants or in targeted regions) and annotates them with the data from dbSnp, Clinvar, phastCons and dbNSFP databases. It creates a temporary directory in the output directory where the files from each step are stored.
# The script takes 5 arguments:
#  v = reference genome version (either hg19 or hg38)
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


# setting the parameters
REFSEQ="/"$version/$version".fa"
SNPEFF="java -jar /snpEff/snpEff.jar"
SNPSIFT="java -jar /snpEff/SnpSift.jar"
DBSNP="/dbSNP/00-All.vcf.gz" # the same file for hg19 and hg38 
CLINVAR="/"$version"/ClinVar/clinvar.vcf.gz"
PHASTCONS="/"$version"/phastCons/"
DBNSFP="/"$version"/dbNSFP/dbNSFP.txt.gz"
GWASCAT="/GWAS_catalog/gwascatalog.txt"
VT="/vt"



echo "Genome version: $version"
echo "Bed file with the regions of interest: $regions"


for dir in $(ls -1 $in_dir | grep -e "serija"); do

  echo processing $dir :
  
  # create the temporary directory
  mkdir $out_dir/$dir/
  mkdir $out_dir/$dir/tmp
  temp=$out_dir/$dir/tmp/
  name=$dir

  ### variant calling
  bams=$(find $in_dir/$dir/ | grep -e ".bam$")
  
  if [ $regions == "false" ] 
  then
    
    # this command calls all variants 
    echo "Calling and indexing variants..."
    freebayes -f $REFSEQ --haplotype-length -1 --min-repeat-entropy 0 $bams | bcftools sort > $temp/$name.sorted.vcf 2> $temp/fb_error_log.txt
    $VT decompose -s $temp/$name.sorted.vcf | $VT normalize -n -r $REFSEQ - | $VT decompose_blocksub - -o $temp/$name.sorted.decomp.norm.vcf 2> $temp/normalisation_error_log.txt
    bgzip $temp/$name.sorted.decomp.norm.vcf
    tabix -p vcf $temp/$name.sorted.decomp.norm.vcf.gz
  
  else
      
    # this command calls variants only in regions of interest/gene panel (provided BED file)
    echo "Calling and indexing variants only in the targeted regions..."
    freebayes -f $REFSEQ -t $regions --haplotype-length -1 --min-repeat-entropy 0 $bams | bcftools sort > $temp/$name.sorted.vcf 2> $temp/fb_error_log.txt
    $VT decompose -s $temp/$name.sorted.vcf | $VT normalize -n -r $REFSEQ - | $VT decompose_blocksub - -o $temp/$name.sorted.decomp.norm.vcf 2> $temp/norm_decomp_error_log.txt
    bgzip $temp/$name.sorted.decomp.norm.vcf
    tabix -p vcf $temp/$name.sorted.decomp.norm.vcf.gz

  fi
  
  
  
  
  ### variant annotation
  # SnpEff: add ANN field 
  echo "annotating variants with the SnpEff..."
  echo processing $temp/$name.sorted.decomp.norm.vcf
  $SNPEFF -v -noStats $version $temp/$name.sorted.decomp.norm.vcf > $temp/$name.sorted.decomp.norm.ann.vcf 2> $temp/snpeff_error_log.txt

  
  #SnpSift: add information from various databases
  echo "SnpSift: annotating variants with the dbSnp, Clinvar and phastCons databases..."
  echo processing $temp/$name.sorted.decomp.norm.ann.vcf
  bcftools sort $temp/$name.sorted.decomp.norm.ann.vcf | $SNPSIFT annotate $DBSNP | $SNPSIFT annotate $CLINVAR | $SNPSIFT phastCons $PHASTCONS - > $temp/$name.sorted.decomp.norm.ann.dbSnp.Clinvar.phastCons.vcf 2> $temp/error_log_snpsift1.txt
  
  # this part is separated because the dbNSFP annotation won't work in the pipeline 
  echo "SnpSift: annotating variants with the dbNSFP scores and adding variant type and zygosity information..."
  echo processing $temp/$name.sorted.decomp.norm.ann.dbSnp.Clinvar.phastCons.vcf
  $SNPSIFT dbnsfp -db $DBNSFP $temp/$name.sorted.decomp.norm.ann.dbSnp.Clinvar.phastCons.vcf | $SNPSIFT varType - | $SNPSIFT gwasCat -db $GWASCAT - > $out_dir/$name.sorted.decomp.norm.ann.dbSnp.Clinvar.phastCons.dbNSFP.gwas.vcf 2> $temp/error_log_snpsift2.txt
  
  
  
  ### deleting temporary files
  if [ $keep == "false" ] 
  then
   
   echo "Deleting temporary files..."
   rm -rf $temp     
   
  fi



done






