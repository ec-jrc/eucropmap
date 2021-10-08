#!/bin/bash
##Go to working dir
#cd /eos/jeodpp/data/projects/REFOCUS/cropclassif

## Argument one passed by the HTCondor submit file
# List file with the path to the imput raster
LIST_product=$1

## Argument two passed by the HTCondor submit file
# Proccess index used for read the line number
N_LINE=$3
##Becouse the HTCondor index starts on 0 we must add one in order to read the line 1
N_LINE=$((N_LINE + 1))


## Argument three passed by the HTCondor submit file
# Similar to the Array index
CLUSTER=$2
echo $LIST_product
echo "Line number: "$N_LINE
echo "Domain: "$UID_DOMAIN
echo "Cluster: "$CLUSTER
echo "Updated"
###OUTPUT FOLDER
OUTPUT_folder=$4
echo "Ouput folder: "$OUTPUT_folder
#Read line number according with the proccess index
INPUT_product=$(awk "NR==$N_LINE" $LIST_product) 
echo "INPUT_product: " $INPUT_product 
FILENAME=$(basename "$INPUT_product")
echo "File name: " $FILENAME
PATHFOLDER=$(dirname "$INPUT_product")
echo "Folder path: " $PATHFOLDER

outdir=/scratch/job${CLUSTER}_${N_LINE}/
mkdir -p $outdir
echo "working dir: " $outdir
ls $outdir



echo python3 s1_classify_2levels.py $INPUT_product $outdir "$5" "$6"
python3 s1_classify_2levels.py $INPUT_product $outdir "$5" "$6"
#echo pip3 freeze > $OUTPUT_folder'pip-freeeze.txt'

mv $outdir* $OUTPUT_folder
rm -rf $outdir
