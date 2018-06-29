#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Define variables with directories
wd=$1
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

OUTPUT_DIR=$2
PROJECT_NAME=$3

if [ ! -d ${OUTPUT_DIR}bsubLogs ]; then
    mkdir ${OUTPUT_DIR}bsubLogs
else
    a=`ls -1 ${OUTPUT_DIR}bsubLogs/ | wc -l`
    if [ $a -ne 0 ]; then
        rm ${OUTPUT_DIR}bsubLogs/*
    fi
fi

if [ ! -d ${OUTPUT_DIR}keep ]; then
     mkdir ${OUTPUT_DIR}keep/
fi

KEEP=${OUTPUT_DIR}keep/
python ${scrd}randomizeFam.py ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.fam $KEEP 500
    
z=0
for k in $(find $KEEP -type f); do
    let z++
    bsub -P hla -q week -o ${OUTPUT_DIR}bsubLogs/hlaImpute.$z.out -R "rusage[mem=10]" "sh ${hlad}SNP2HLA/SNP2HLA.sh ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS ${hlad}T1DGC/T1DGC_REF_exome ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.tmp.${z}.HLA $k"
done

check=0
while [ $check -eq 0 ]; do
    q=`ls -1 ${OUTPUT_DIR}bsubLogs/ | grep "hlaImpute" | wc -l | awk '{print $1}'`
    if [ $(($q)) -eq $z ]; then
        check=1
        break
    else
        date=`date +"%m/%d/%Y %T"`                    
        sleep 10m
    fi
done

FILES=${OUTPUT_DIR}bsubLogs/
exit=0
for f in $(find $FILES -type f); do
    p=`grep -c "Successfully completed." $f`
    if [ $(($p)) -eq 0 ]; then
        echo "Log file ${f} failed... exiting pipeline"
        exit=1
    fi
done

if [ $exit -eq 1 ]; then
    exit 1
fi

for i in $( seq 2 $z ); do
            echo "${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.tmp.${i}.HLA.bed ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.tmp.${i}.HLA.bim ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.tmp.${i}.HLA.fam" >> ${OUTPUT_DIR}mergeListHLA.txt
done


plink --bfile ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.tmp.1.HLA --merge-list ${OUTPUT_DIR}mergeListHLA.txt --make-bed --out ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.HLA
python ${scrd}mergeDosage.py ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.tmp $z > ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.HLA.dosage

echo "" > ${OUTPUT_DIR}hla.checksum