#!/bin/bash

for i in {1..5}
do

  run=run$i
  mkdir /dir/$run
  
  echo "currently processing" $run
  
  bash /script.sh -v hg38 -i /.../bams/$run -r "/.../amplicons.bed" -k true -F 0.01 -o /dir/$run/ 
  
done
