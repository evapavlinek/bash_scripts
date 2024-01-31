#!/bin/bash

#Perform the analysis of depth and breadth of coverage for target regions using bed tools

files=($(find ./OUTPUT/ -name '*.sorted.bam'))

#calculate the coverage
parallel --jobs 8 --ungroup 'echo processing {}; fname=$(basename {}); bedtools genomecov -ibam {} -g /hgXX.genome -dz > /dir/$fname.cov' ::: ${files[@]} 


