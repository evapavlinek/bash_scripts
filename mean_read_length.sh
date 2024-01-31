#!/bin/bash

echo "collecting alignment summary metrics for the bam files..."
bams=($(find /bams/ -name '*.sorted.filtered.bam'))
parallel --jobs 6 --ungroup 'echo processing {}; java -jar /picard.jar CollectAlignmentSummaryMetrics -R /hgXX.fa -I {} -O "/dir/"$(basename {} ".sorted.filtered.bam")_alignment_summary_metrics.txt' ::: ${bams[@]} 2> error_log.txt

