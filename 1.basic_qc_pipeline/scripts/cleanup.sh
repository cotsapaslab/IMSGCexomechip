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

## Check for Errors
function die(){
    echo "Error in cleanup.sh..."
    exit 1
    }
trap die ERR

CLEANUP=`grep "CLEANUP:" $META | awk '{print $2}'`
if [ $CLEANUP -eq 1 ]; then
    date=`date +"%m/%d/%Y %T"`
    echo "[$i] Cleanup -- $date" >> $LOG; let i++

    if [ ! -d ${OUTPUT_DIR}logs ]; then
        mkdir ${OUTPUT_DIR}logs/
        mv ${OUTPUT_DIR}*.log ${OUTPUT_DIR}logs
    fi

    if [ ! -d ${OUTPUT_DIR}keep ]; then
        mkdir ${OUTPUT_DIR}keep/
    fi

    if [ ! -d ${OUTPUT_DIR}plots ]; then
        mkdir ${OUTPUT_DIR}plots/
        mv ${OUTPUT_DIR}*.png ${OUTPUT_DIR}plots/
    fi

    rm ${OUTPUT_DIR}*.SSD*
    rm ${OUTPUT_DIR}*.Info*

    rm ${OUTPUT_DIR}*.nosex
    rm ${OUTPUT_DIR}*.ped
    rm ${OUTPUT_DIR}*.map
    rm ${OUTPUT_DIR}*.hh
    rm ${OUTPUT_DIR}*.geno
    rm ${OUTPUT_DIR}*.snp
    rm ${OUTPUT_DIR}*.ind
    rm ${OUTPUT_DIR}*.pedind            
    rm ${OUTPUT_DIR}tmp.*
    rm ${OUTPUT_DIR}par*
    rm ${OUTPUT_DIR}z*
    rm ${OUTPUT_DIR}gencall*    
    rm ${OUTPUT_DIR}y*
    rm ${OUTPUT_DIR}merge*
    rm ${OUTPUT_DIR}MERGE*
    gzip ${OUTPUT_DIR}genome.*genome
    mv ${OUTPUT_DIR}keep* ${OUTPUT_DIR}keep/
    mv ${OUTPUT_DIR}offspring* ${OUTPUT_DIR}keep/
    mv ${OUTPUT_DIR}dups* ${OUTPUT_DIR}keep/
    mv ${OUTPUT_DIR}sibs* ${OUTPUT_DIR}keep/
    mv ${OUTPUT_DIR}piHat* ${OUTPUT_DIR}keep/
    mv ${OUTPUT_DIR}bad* ${OUTPUT_DIR}keep/        
    tar -pczf ${OUTPUT_DIR}keptSNPsSamples.tar.gz ${OUTPUT_DIR}keep/
fi
