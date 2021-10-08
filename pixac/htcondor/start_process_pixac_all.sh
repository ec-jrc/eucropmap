#!/bin/bash
##Go to working dir
#cd /eos/jeodpp/data/projects/REFOCUS/cropclassif

## Argument one passed by the HTCondor submit file
# List file with the path to the input raster
LIST_product=${3}
echo "input file":$LIST_product

## Argument two passed by the HTCondor submit file
# Proccess index used for read the line number
N_LINE=${2}
##Becouse the HTCondor index starts on 0 we must add one in order to read the line 1
N_LINE=$(($N_LINE + 1))


## Argument three passed by the HTCondor submit file
# Similar to the Array index
CLUSTER=${1}

echo "Line number: "$N_LINE
#echo "Domain: "$UID_DOMAIN
echo "Cluster: "$CLUSTER
#echo "Updated"
echo "INPUT_product: " $LIST_product
FILENAME=$(basename "$LIST_product")
echo "File name: " $FILENAME
PATHFOLDER=$(dirname "$LIST_product")
echo "Folder path: " $PATHFOLDER

###OUTPUT FOLDER
OUTPUT_folder=${4}
echo "Ouput folder: "$OUTPUT_folder

classif_dir=${5}
echo "classif_dir: "$classif_dir

modelpath=${6}
echo "modelpath: "$modelpath

csv_file=${7}
echo "csv_file: "$csv_file

stratum=${8}
echo "stratum: "$stratum

outdir=/scratch/REFOCUS/job${CLUSTER}_${N_LINE}/
mkdir -p $outdir
echo "working dir: " $outdir
ls $outdir

echo "python3 /eos/jeodpp/data/projects/REFOCUS/data/S1_GS/v7/pixaccuracy/pixac/run_eucropmap_all_htc.py $LIST_product $outdir $classif_dir $modelpath $csv_file $stratum"
python3 /eos/jeodpp/data/projects/REFOCUS/data/S1_GS/v7/pixaccuracy/pixac/run_eucropmap_all_htc.py $LIST_product $outdir $classif_dir $modelpath $csv_file $stratum

#echo "python3 /eos/jeodpp/data/projects/REFOCUS/data/S1_GS/v7/pixaccuracy/pixac/run_eucropmap_all_htc.py $INPUT_product $outdir $classif_dir $modelpath $csv_file" #$stratum"
#python3 /eos/jeodpp/data/projects/REFOCUS/data/S1_GS/v7/pixaccuracy/pixac/run_eucropmap_all_htc.py $INPUT_product $outdir $classif_dir $modelpath $csv_file #$stratum
#echo pip3 freeze > $OUTPUT_folder'pip-freeeze.txt'

mv $outdir* $OUTPUT_folder
rm -rf $outdir

#sleep infinity
echo "done"