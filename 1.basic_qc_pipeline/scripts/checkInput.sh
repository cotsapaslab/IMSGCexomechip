#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Define variables
META=$1
i=$2
LOG=$3

OUTPUT_DIR=`grep "Output Directory:" $META | awk '{print $3}'`
PHENO=`grep "Phenotype File:" $META | awk '{print $3}'`
COHORT=`grep "Cohort File:" $META | awk '{print $3}'`
GAP_SUMMARY=`grep "GAP Summary File:" $META | awk '{print $4}'`
PROJECT_NAME=`grep "Project Name:" $META | awk '{print $3}'`

GENCALL=`grep "gencall file:" $META | awk '{print $3}'`
ZCALL=`grep "zcall file:" $META | awk '{print $3}'`

## Define variables with directories
wd=`grep "Script Directory:" $META | awk '{print $3}'`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

## Check for Errors
function die(){
    echo "Error in checkInput.sh..."
    exit 1
    }
trap die ERR

date=`date +"%m/%d/%Y %T"`
echo "[$i] Check Whether Inputs Exist -- $date" >> $LOG; let i++

## Check for output directory existence
if [ ! -d $OUTPUT_DIR ]; then
    echo "Output directory $OUTPUT_DIR does not exist." >> $LOG; exit=1
fi

## Check for plink files existence
exists=0    
for root in $GENCALL $ZCALL; do
    if [ -n "$root" ]; then
        for ext in bed bim fam; do
            if [ ! -e $root.$ext ]; then
                echo "File does not exist: $root.$ext" >> $LOG; exit=1
            else
                exists=1
            fi            
        done
    fi
done

    ## If no existing inputs, exit
    if [ $exists -eq 0 ]; then
        echo "No Input Files exist; Make sure meta file is formatted correctly" >> $LOG
        exit 1
    fi

    ## Check for plink file and output directory existence errors
    if [ "$exit" = "1" ]; then
        date=`date +"%m/%d/%Y %T"`
        echo "" >> $LOG
        echo "Exiting pipeline ... $date" >> $LOG    
        exit 1
    fi

## Check for _C appended to HapMap Controls
for root in $GENCALL $ZCALL; do
    if [ -n "$root" ]; then
        hapmap=`python ${scrd}checkForHapMapAliquot.py $root.fam | wc -l | awk '{ print $1 }'`        
        if [ $hapmap -ne 0 ]; then
            echo "HapMap samples have aliquots (ex: _C) appended  in ${root}.fam -- Please remove from IID and rerun scripts." >> $LOG; exit=1
            python ${scrd}checkForHapMapAliquot.py $root.fam >> $LOG
            echo "" >> $LOG
        fi
    fi
done

## Check if any failures
if [ "$exit" = "1" ]; then
    date=`date +"%m/%d/%Y %T"`    
    echo "Exiting pipeline ... $date" >> $LOG
    exit 1
fi

## Check for dups
for root in $GENCALL $ZCALL; do
    if [ -n $root ]; then
        awk '{print $2 }' $root.fam | sort -gk 1 | uniq -c | awk '{ if ($1 > 1) print $2}' > ${OUTPUT_DIR}dups.txt
        
        nDups=`wc -l ${OUTPUT_DIR}dups.txt | awk '{print $1}'`
        if [ $nDups -ne 0 ]; then
            echo "Duplicates present in ${root}.fam - Please remove using [ ${scrd}removeDuplicates.sh ] and restart pipeline" >> $LOG
            cat ${OUTPUT_DIR}dups.txt >> $LOG
            echo "" >> $LOG
            exit=1
        fi
        rm ${OUTPUT_DIR}dups.txt            
    fi
done

## Check if any failures
if [ "$exit" = "1" ]; then
    date=`date +"%m/%d/%Y %T"`    
    echo "Exiting pipeline ... $date" >> $LOG
    rm ${OUTPUT_DIR}tmp.fam
    exit 1
fi

## Check alleles are in FHG19 -- if not convert them
for root in $GENCALL $ZCALL; do
        
    ## Update alleles from AB to ATGC        
    ab=`awk '{ if ($5 == "B" || $6 == "B") print $0 }' $root.bim | wc -l | awk '{ print $1 }'`
    if [ "$(($ab))" -ne "0" ]; then
        echo "Alleles for $root are in AB format. Use the allele map in [ ${supd}Allele_Map_AB_to_HG19_ID.txt ] to convert to FHG19 and restart pipeline." >> $LOG
        exit=1
    fi
            
    ## Update alleles from ID to ATGC,- 
    id=`awk '{ if ($5 == "D" || $6 == "D" || $5 == "I" || $6 == "I") print $0 }' $root.bim | wc -l | awk '{ print $1 }'`
    if [ "$(($id))" -ne "0" ]; then
        echo "Indels for $root are in I/D format. Use the allele map in [ ${supd}Allele_Map_ID_to_FHG19.txt ] to convert." >> $LOG
        exit=1
    fi
            
    ## Update alleles from HG19 to FHG19 ##Need to test later 1/4/13
    flip=`python ${scrd}checkIfHG19.py $root.bim ${supd}fhg19.bim | wc -l | awk '{ print $0 }'`
    if [ "$(($flip))" -ne "0" ]; then
        echo "Alleles are in HG19 format. Use the allele map in [ ${supd}Allele_Map_HG19_to_FHG19.txt ] and the flip list in [ ${supd}flipBOT.txt ] to convert." >> $LOG
        exit=1
    fi
done

## Check if any failures
if [ "$exit" = "1" ]; then
    date=`date +"%m/%d/%Y %T"`    
    echo "Exiting pipeline ... $date" >> $LOG
    exit 1
else
    date=`date +"%m/%d/%Y %T"`
    echo "All required inputs exist and are in correct format ... $date" >> $LOG
    echo "" >> $LOG
    echo "" >> $LOG
fi        