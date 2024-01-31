#!/bin/bash

# this command calls all variants 
echo "calling variants..."
bams=($(find /dir/ -name '*.sorted.filtered.marked_dup.bam'))
parallel --jobs 6 --ungroup 'echo processing $(basename {}); freebayes -f /hgXX.fa {} > "/dir/"$(basename {} .bam).vcf' ::: ${bams[@]} 2> error_log.txt

# this command calls variants only in regions of interest/gene panel (provided BED file)
echo "calling variants only in the targeted regions..."
parallel --jobs 6 --ungroup 'echo processing $(basename {}); freebayes -f /hgXX.fa -t /dir/probes.bed {} > "/dir/"$(basename {} .bam).targeted_regions.vcf' ::: ${bams[@]} 2> error_log.txt

