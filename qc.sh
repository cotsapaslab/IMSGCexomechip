#!/bin/bash -E
set -e
alias plink='plink --noweb --silent --allow-no-sex'

## Define variables from Meta File
META=$1

OUTPUT_DIR=`grep "Output Directory:" $META | awk '{print $3}'`
PHENO=`grep "Phenotype File:" $META | awk '{print $3}'`
COHORT=`grep "Cohort File:" $META | awk '{print $3}'`
GAP_SUMMARY=`grep "GAP Summary File:" $META | awk '{print $4}'`
PROJECT_NAME=`grep "Project Name:" $META | awk '{print $3}'`

## Define directories
wd=`grep "Script Directory:" $META | awk '{print $3}'`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

exit=0
if [ ! -d ${scrd} ]; then
    echo "${scrd} does not exist"; exit=1   
fi
if [ ! -d ${supd} ]; then
    echo "${supd} does not exist"; exit=1
fi
if [ ! -d ${eigd} ]; then
    echo "${eigd} does not exist"; exit=1
fi
if [ ! -d ${hlad} ]; then
    echo "${hlad} does not exist"; exit=1
fi
if [ ! -d ${cond} ]; then
    echo "${cond} does not exist"; exit=1
fi

if [ $exit -eq 1 ]; then
    exit 1
fi

## Initialize Log File
LOG=`sh ${scrd}makeLogFile.sh $META`
i=1

## Check for Errors
function die(){
    echo "Exited pipeline." >> $LOG
    echo "Exited pipeline. See $LOG for more details."        
    exit 1
    }
trap die ERR

exit=0

## Execute shell scripts based on options in meta file
CHECK_INPUT=1 # Set to 1 by default
if [ $CHECK_INPUT -eq 1 ]; then
    sh ${scrd}checkInput.sh $META $i $LOG
    let i++
fi

SAMPLE_QC=`grep "SAMPLE_QC:" $META | awk '{print $2}'`
if [ $SAMPLE_QC -eq 1 ]; then
    sh ${scrd}sampleQC.sh $META $i $LOG
    let i++
fi

RELATEDNESS_CHECK=`grep "RELATEDNESS_CHECK:" $META | awk '{print $2}'`
if [ $RELATEDNESS_CHECK -eq 1 ]; then
    sh ${scrd}relatednessCheck.sh $META $i $LOG
    let i++
fi

SEX_CHECK=`grep "SEX_CHECK:" $META | awk '{print $2}'`
if [ $SEX_CHECK -eq 1 ]; then
    sh ${scrd}sexCheck.sh $META $i $LOG
    let i++
fi

MENDEL_CHECK=`grep "MENDEL_CHECK:" $META | awk '{print $2}'`
if [ $MENDEL_CHECK -eq 1 ]; then
    sh ${scrd}mendelCheck.sh $META $i $LOG
    let i++
fi
            
PCA_ROUND1=`grep "PCA_ROUND1:" $META | awk '{print $2}'`
if [ $PCA_ROUND1 -eq 1 ]; then
    sh ${scrd}pcaRound1.sh $META $i $LOG
    let i++
fi

SITE_QC=`grep "SITE_QC:" $META | awk '{print $2}'`
if [ $SITE_QC -eq 1 ]; then
    sh ${scrd}siteQC.sh $META $i $LOG
    let i++
fi

date=`date +"%m/%d/%Y %T"`
echo "Finished pipeline ... $date" >> $LOG
echo "---------------------------------------------------------" >> $LOG