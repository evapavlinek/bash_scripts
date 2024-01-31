#!/bin/bash

#Perform the analysis of depth and breadth of coverage for target regions using bed tools

files=($(find ./OUTPUT/ -name '*.sorted.bam'))

#calculate the coverage
parallel --jobs 8 --ungroup 'echo processing {}; fname=$(basename {}); bedtools coverage -a ./probe_regions.bed -b {} > /dir/$fname.cov' ::: ${files[@]} 


#calculate the coverage with the low mapping quality reads removed
parallel --jobs 8 --ungroup 'echo processing {}; fname=$(basename {}); samtools view -b -q 10 {} > /dir/filtered/$fname.filter' ::: ${files[@]} 


filtered=($(find ./dir/filtered/ -name '*.filter'))

parallel --jobs 8 --ungroup 'echo processing {}; fname=$(basename {}); bedtools coverage -a ./probe_regions.bed -b {} > /dir/filtered/$fname.cov' ::: ${filtered[@]} 

