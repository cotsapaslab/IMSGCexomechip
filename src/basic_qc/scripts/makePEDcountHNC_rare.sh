#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Variables from stdin
OUTPUT_DIR=$1
root=$2
z=$3
rareSNPs=$4
scrd=$5
i=$6
kg=$7

if [ $kg -eq 1 ]; then
    plink --bfile $root --extract $rareSNPs --recode --out ${OUTPUT_DIR}tmp.rare.${z} --keep ${OUTPUT_DIR}tmp.list`printf "%03i\n" $i`
else
    plink --bfile $root --geno 0.05 --max-maf 0.05 --recode --out ${OUTPUT_DIR}tmp.rare.${z} --keep ${OUTPUT_DIR}tmp.list`printf "%03i\n" $i`
fi

python ${scrd}countHetNoCalls.py ${OUTPUT_DIR}tmp.rare.${z}.ped > ${OUTPUT_DIR}tmp.rare.${z}.hetnocall
