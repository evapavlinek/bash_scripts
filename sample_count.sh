#!/bin/bash

for dir in /dir/*; do
    
  dir_name=$(basename $dir)
  input_dir=/dir/"$dir_name"/
  
  echo $dir_name
  
  cd $dir
  echo $(ls | wc -l)
  
  cd $input_dir
  echo $(ls | wc -l)
  
  echo ""

done

