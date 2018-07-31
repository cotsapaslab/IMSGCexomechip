#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Define variables
OUTPUT_DIR=$1
PROJECT_NAME=$2
SCRIPT_DIR=$3

## Define variables with directories
wd=$SCRIPT_DIR
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

array=( $@ )

if [ -e ${OUTPUT_DIR}mergeList.txt ]; then
    rm ${OUTPUT_DIR}mergeList.txt
fi

for i in $(seq 3 $((${#array[@]} - 1))); do
    echo ${array[$i]}.bed ${array[$i]}.bim ${array[$i]}.fam >> ${OUTPUT_DIR}mergeList.txt
done

plink --bfile ${cond}controls --merge-list ${OUTPUT_DIR}mergeList.txt --make-bed --out ${OUTPUT_DIR}RAW_MERGE

python ${scrd}findIntersectionSNPs.py ${cond}controls.bim ${array[@]:3:${#array[@]}} > ${OUTPUT_DIR}keepMergeSNPs.txt

plink --bfile ${OUTPUT_DIR}RAW_MERGE --extract ${OUTPUT_DIR}keepMergeSNPs.txt --keep ${OUTPUT_DIR}european.txt --hwe 1e-5 --make-bed --out ${OUTPUT_DIR}MERGE.eur

awk '{ print $2 }' ${OUTPUT_DIR}MERGE.eur.bim > ${OUTPUT_DIR}keepMergeSNPs2.txt
plink --bfile ${OUTPUT_DIR}RAW_MERGE --extract ${OUTPUT_DIR}keepMergeSNPs2.txt --make-bed --out ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS --geno 0.005

if [ ! -e ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.fam ]; then
    echo "Problem merging... check alleles match."
    exit 1
fi

nSamples=`wc -l ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.fam | awk '{print $1 }'`
nSNPS=`wc -l ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.bim | awk '{print $1 }'`
echo "After merging all batches with the controls set, there are $nSamples samples and $nSNPS snps." 
for j in 1 2; do
    echo "" 
done

