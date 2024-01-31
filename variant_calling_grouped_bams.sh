#!/bin/bash

echo "calling variants in the grouped bam files"
bams=($(find "/dir/" -name '*.sorted.filtered.rg.bam'))
freebayes -f /hgXX.fa $bams > "./grouped_bams.vcf"
