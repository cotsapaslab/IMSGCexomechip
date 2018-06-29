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

GC_CR=`grep "gencall_missrate:" $META | awk '{print $2}'`
DROP_SEX=`grep "removesexerror:" $META | awk '{print $2}'`

date=`date +"%m/%d/%Y %T"`
echo "[$i] Sex Check -- $date" >> $LOG; ((i++))
plink --bfile $GENCALL --check-sex --out ${OUTPUT_DIR}sexCheck --geno $GC_CR
awk '{ if ($3 == -9) print $1, $2, $4 }' ${OUTPUT_DIR}sexCheck.sexcheck > ${OUTPUT_DIR}sexUpdate.txt
awk '{ if ($3 != -9 && $4 != 0 && $5 == "PROBLEM" && $3 != "N/A" && $3 != 0) print $0 }' ${OUTPUT_DIR}sexCheck.sexcheck > ${OUTPUT_DIR}sexCheckErrors.txt
nError=`wc -l ${OUTPUT_DIR}sexCheckErrors.txt | awk '{print $1 }'`
if [ $(($nError)) -ne 0 ]; then
        echo "There are $nError sex check errors." >> $LOG
        cat ${OUTPUT_DIR}sexCheckErrors.txt >> $LOG
        if [ $DROP_SEX = "1" ]; then
            echo "These samples were dropped from final output." >> $LOG
            awk '{ print $1, $2 }' ${OUTPUT_DIR}sexCheckErrors.txt >> ${OUTPUT_DIR}bad.samples.txt
        fi
fi

for j in 1 2; do
    echo "" >> $LOG
done    

echo "$date" > ${OUTPUT_DIR}sexcheck.checksum
