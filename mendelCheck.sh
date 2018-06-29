#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Define variables
META=$1
i=$2
LOG=$3

GENCALL=`grep "gencall file:" $META | awk '{print $3}'`
ZCALL=`grep "zcall file:" $META | awk '{print $3}'`

OUTPUT_DIR=`grep "Output Directory:" $META | awk '{print $3}'`
PHENO=`grep "Phenotype File:" $META | awk '{print $3}'`
COHORT=`grep "Cohort File:" $META | awk '{print $3}'`
GAP_SUMMARY=`grep "GAP Summary File:" $META | awk '{print $4}'`
PROJECT_NAME=`grep "Project Name:" $META | awk '{print $3}'`

## Define variables with directories
wd=`grep "Script Directory:" $META | awk '{print $3}'`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

## Check for Errors
function die(){
    echo "Error in mendelCheck.sh..."
    exit 1
    }
#trap die ERR


GC_CR=`grep "gencall_missrate:" $META | awk '{print $2}'`
DROP_MENDEL=`grep "removemendelerror:" $META | awk '{print $2}'`
nmendel=`grep "maxMendelError:" $META | awk '{print $2 }'` 

date=`date +"%m/%d/%Y %T"`
echo "[$i] Mendel Check -- $date" >> $LOG; ((i++))

exists=`awk '{ if ($3 != 0 || $4 != 0) print $0 }' $GENCALL.fam | wc -l | awk '{ print $1 }'`
if [ $(($exists)) -ne 0 ]; then
    plink --bfile $GENCALL --mendel --out ${OUTPUT_DIR}mendelCheck --geno $GC_CR --maf 0.05
    if [ $DROP_MENDEL = "1" ]; then
        awk -v nmendel=$nmendel '{ if ($5 > nmendel) print $1}' ${OUTPUT_DIR}mendelCheck.fmendel > ${OUTPUT_DIR}familyIDs.todrop.txt    
        python ${scrd}dropMendelErrors.py $GENCALL.fam ${OUTPUT_DIR}familyIDs.todrop.txt > ${OUTPUT_DIR}mendelErrors.txt
        cat ${OUTPUT_DIR}mendelErrors.txt >> ${OUTPUT_DIR}bad.samples.txt
        echo "Dropped mendel errors. See ${OUTPUT_DIR}mendelErrors.txt for samples dropped" >> $LOG
    fi
else
    echo "No families defined in [ $GENCALL.fam ]. Stopping Mendel Check..." >> $LOG
fi

for j in 1 2; do
    echo "" >> $LOG
done    

echo "$date" > ${OUTPUT_DIR}mendelcheck.checksum
