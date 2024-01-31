#!/bin/bash

#  Direcotry with the annotated files
INPUT_DIR=/dir

# Filter to be used for selecting files from 
FILTER="filters.RDS"

#get all the files that end in filter from input directory
arr=($INPUT_DIR/*$FILTER)

# iterate through array using a counter
for ((i=0; i<${#arr[@]}; i++)); do
	    #do something to each element of array
	    SAMPLE="${arr[$i]}"
		SAMPLE=$(basename "$SAMPLE")
		SAMPLE=${SAMPLE%"$FILTER"}
		SCRIPT="Rscript -e 'library(rmarkdown); rmarkdown::render("Rscript_name.R ${SAMPLE}")' >${SAMPLE}_log.txt &"
		#now to run the script 
		echo "knitting pdf for sample $SAMPLE" 
                $SCRIPT "$SAMPLE"
		
	done
