#!/bin/bash

for dir in /dir/*; do
  
  if [ -z "$(ls -A $dir)" ]; then
    echo "$dir is empty"
  else
    echo "$dir is not Empty"
  fi

done


