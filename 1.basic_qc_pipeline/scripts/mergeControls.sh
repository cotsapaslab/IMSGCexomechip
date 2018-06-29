#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Define variables
META=$1
i=$2
LOG=$3

GENCALL1=`grep "GenCall PLINK file1 root path:" $META | awk '{print $6}'`
ZCALL1=`grep "zCall PLINK file1 root path:" $META | awk '{print $6}'`
GENCALL2=`grep "GenCall PLINK file2 root path:" $META | awk '{print $6}'`
ZCALL2=`grep "zCall PLINK file2 root path:" $META | awk '{print $6}'`
GENCALL3=`grep "GenCall PLINK file3 root path:" $META | awk '{print $6}'`
ZCALL3=`grep "zCall PLINK file3 root path:" $META | awk '{print $6}'`
QC_IN=`grep "qc'd dataset root path:" $META | awk '{print $5}'`

OUTPUT_DIR=`grep "Output Directory:" $META | awk '{print $3}'`
PHENO=`grep "Phenotype File:" $META | awk '{print $3}'`
COHORT=`grep "Cohort File:" $META | awk '{print $3}'`
GAP_SUMMARY=`grep "GAP Summary File:" $META | awk '{print $4}'`
PROJECT_NAME=`grep "Project Name:" $META | awk '{print $3}'`

## Define variables with directories
wd=`pwd`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

MERGE_CONTROLS=`grep "MERGE_CONTROLS:" $META | awk '{print $2}'`
if [ $MERGE_CONTROLS -eq 1 ]; then
    date=`date +"%m/%d/%Y %T"`
    echo "[$i] Merge with Combined Controls -- $date" >> $LOG; ((i++))

    n=0
    controls=${cond}controls
    if [ -n $QC_IN ];
        IN=$QC_IN
    else
        IN=${OUTPUT_DIR}${PROJECT_NAME}
    fi
    
    for root in $IN; do
        ((n++))
        plink --bfile $IN --bmerge $controls.bed $controls.bim $controls.fam --make-bed --out ${OUTPUT_DIR}MERGE.$n
        MERGE=${OUTPUT_DIR}MERGE.$n
    done

    python ${scrd}findIntersectionSNPs.py $controls.bim $IN.bim > ${OUTPUT_DIR}keepMergeSNPs.txt

    cat ${OUTPUT_DIR}bad.samples.*.txt ${OUTPUT_DIR}bad.samples.sex_mendel.txt > ${OUTPUT_DIR}finalBadSamples.txt
    cat ${OUTPUT_DIR}european.*.txt ${cond}european.txt > ${OUTPUT_DIR}keepEuropeanSamples.txt
    plink --bfile $MERGE --extract ${OUTPUT_DIR}keepMergeSNPs.txt --keep ${OUTPUT_DIR}keepEuropeanSamples.txt --hwe 1e-5 --make-bed --out ${OUTPUT_DIR}MERGE.eur --remove ${OUTPUT_DIR}finalBadSamples.txt
    awk '{ print $2 }' ${OUTPUT_DIR}MERGE.eur.bim > ${OUTPUT_DIR}keepMergeSNPs2.txt
    plink --bfile $MERGE --extract ${OUTPUT_DIR}keepMergeSNPs2.txt --make-bed --out ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS --geno 0.005 --remove ${OUTPUT_DIR}bad.samples.sex_mendel.txt

    if [ ! -e ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.fam ]; then
        echo "Problem merging... check alleles match." >> $LOG
        exit 1
    fi

    nSamples=`wc -l ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.fam | awk '{print $1 }'`
    nSNPS=`wc -l ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.bim | awk '{print $1 }'`
    echo "After merging with the controls set, there are $nSamples samples and $nSNPS snps." >> $LOG
    for j in 1 2; do
        echo "" >> $LOG
    done
fi
