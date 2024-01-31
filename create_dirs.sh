#!/bin/bash

for dir in $(ls -1 /dir/ | grep -e "pattern"); do
  mkdir /dir/"$dir"/
done


