#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Variables from stdin
META=$1
i=$2
LOG=$3

minHet=`grep "minHet:" $META | awk '{print $2}'`
maxHet=`grep "maxHet:" $META | awk '{print $2}'`
maxNoCall=`grep "maxNoCall:" $META | awk '{print $2}'`
use1KGdef=`grep "use1KGdef:" $META | awk '{print $2}'`
use1KGdef=$(($use1KGdef))


# Variables from Meta File
GENCALL=`grep "gencall file:" $META | awk '{print $3}'`
ZCALL=`grep "zcall file:" $META | awk '{print $3}'`

OUTPUT_DIR=`grep "Output Directory:" $META | awk '{print $3}'`
PHENO=`grep "Phenotype File:" $META | awk '{print $3}'`
COHORT=`grep "Cohort File:" $META | awk '{print $3}'`
GAP_SUMMARY=`grep "GAP Summary File:" $META | awk '{print $4}'`
PROJECT_NAME=`grep "Project Name:" $META | awk '{print $3}'`
maxJobs=`grep "maxJobs:" $META | awk '{ print $2 }'`


## Define variables with directories
wd=`grep "Script Directory:" $META | awk '{print $3}'`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/


date=`date +"%m/%d/%Y %T"`
echo "[$i] Begin Sample QC -- $date" >> $LOG; ((i++))
echo "" >> $LOG

if [ ! -d ${OUTPUT_DIR}bsubLogs ]; then
    mkdir ${OUTPUT_DIR}bsubLogs
else
    a=`ls -1 ${OUTPUT_DIR}bsubLogs/ | wc -l | awk '{print $1}'`
    if [ $(($a)) -ne 0 ]; then
        rm ${OUTPUT_DIR}bsubLogs/*
    fi
fi


echo "Splitting up count het calls and no calls..." >> $LOG
gawk '{print $1,$2}' $GENCALL.fam | split -d -a 3 -l 500 - ${OUTPUT_DIR}tmp.list

i=0
a=`ls -1 ${OUTPUT_DIR}tmp.list* | wc -l | awk '{print $1 }'`
a=$(($a - 1))
z=0
commonSNPs=${supd}commonSNPs.txt # Defined from QC'd 1kg snps
rareSNPs=${supd}rareSNPs.txt # Defined from QC'd 1kg snps            
while [ $i -le $a ]; do
    nJobs=`bjobs | wc -l`
    if [ $(($nJobs)) -gt $(($maxJobs)) ]; then
        sleep 5m
        continue
    fi
    
    if [ $(($use1KGdef)) -eq 1 ]; then
        ((z++))
        bsub -P ${PROJECT_NAME} -q hour -W 4:00 -o ${OUTPUT_DIR}bsubLogs/$z.log sh ${scrd}makePEDcountHNC_common.sh $OUTPUT_DIR $GENCALL $z $commonSNPs $scrd $i 1
        ((z++))                    
        bsub -P ${PROJECT_NAME} -q hour -W 4:00 -o ${OUTPUT_DIR}bsubLogs/$z.log sh ${scrd}makePEDcountHNC_rare.sh $OUTPUT_DIR $GENCALL $z $rareSNPs $scrd $i 1

    else
        ((z++))                                    
        bsub -P ${PROJECT_NAME} -q hour -W 4:00 -o ${OUTPUT_DIR}bsubLogs/$z.log sh ${scrd}makePEDcountHNC_common.sh $OUTPUT_DIR $GENCALL $z $commonSNPs $scrd $i 0
        ((z++))                                        
        bsub -P ${PROJECT_NAME} -q hour -W 4:00 -o ${OUTPUT_DIR}bsubLogs/$z.log sh ${scrd}makePEDcountHNC_rare.sh $OUTPUT_DIR $GENCALL $z $rareSNPs $scrd $i 0
    fi
    ((i++))
done

check=0
while [ $check -eq 0 ]; do
    q=`ls -1 ${OUTPUT_DIR}bsubLogs/ | wc -l | awk '{print $1}'`
    if [ $(($q)) -eq $z ]; then
        check=1
        break
    else
        date=`date +"%m/%d/%Y %T"`                    
        echo "Still waiting for jobs to finish... $date" >> $LOG
        sleep 5m
    fi
done

FILES=${OUTPUT_DIR}bsubLogs/
exit=0
for f in $(find $FILES -type f); do
    p=`grep -c "Successfully completed." $f`
    if [ $(($p)) -eq 0 ]; then
        echo "Log file ${f} failed... exiting pipeline" >> $LOG
        echo "Log file ${f} failed... exiting pipeline"
        exit=1
    fi
done

if [ $exit -eq 1 ]; then
    exit 1
fi

date=`date +"%m/%d/%Y %T"`                    
echo "Jobs finished running... $date" >> $LOG

cat ${OUTPUT_DIR}tmp.common.*.hetnocall | awk '{ if (NR == 1 || (!(match($0,"FID")))) print $0 }' > ${OUTPUT_DIR}tmp.common.hetnocall
cat ${OUTPUT_DIR}tmp.rare.*.hetnocall | awk '{ if (NR == 1 || (!(match($0,"FID")))) print $0 }' > ${OUTPUT_DIR}tmp.rare.hetnocall
            
rm ${OUTPUT_DIR}tmp.*.ped
rm ${OUTPUT_DIR}tmp.*.map

## If pheno and cohort file exist, add to counts from above for future plotting color coding
COHORT=`grep "Cohort File:" $META | awk '{print $3}'`
if [ -z "$COHORT" ]; then            
    echo "FID IID cohort" > ${OUTPUT_DIR}cohort.txt
    awk -v pn="$PROJECT_NAME" '{ print $1, $2, pn }' $GENCALL.fam >> ${OUTPUT_DIR}cohort.txt
    COHORT=${OUTPUT_DIR}cohort.txt
fi

PHENO=`grep "Phenotype File:" $META | awk '{print $3}'`
if [ -z "$PHENO" ]; then
    echo "FID IID pheno1" > ${OUTPUT_DIR}pheno.txt
    awk '{ print $1, $2, "-9" }' $GENCALL.fam >> ${OUTPUT_DIR}pheno.txt
    PHENO=${OUTPUT_DIR}pheno.txt
fi
                                
python ${scrd}addPhenotypes.py ${OUTPUT_DIR}tmp.common.hetnocall $COHORT $PHENO > ${OUTPUT_DIR}common.hetnocall
python ${scrd}addPhenotypes.py ${OUTPUT_DIR}tmp.rare.hetnocall $COHORT $PHENO > ${OUTPUT_DIR}rare.hetnocall
            
rm ${OUTPUT_DIR}tmp.*.hetnocall

Rscript ${scrd}plotHetNoCall.r ${OUTPUT_DIR}common.hetnocall ${OUTPUT_DIR}rare.hetnocall ${OUTPUT_DIR}hetnocall.png $minHet $maxHet $maxNoCall

nSamples=`wc -l ${OUTPUT_DIR}rare.hetnocall | awk '{print $1}'`
nSamples=$(($nSamples - 1))

echo "There are $nSamples samples to begin with..." >> $LOG

awk '{ if (($4 < '$minHet' || $4 > '$maxHet') && NR != 1) print $1, $2 }' ${OUTPUT_DIR}rare.hetnocall > ${OUTPUT_DIR}bad.rare.txt
awk '{ if ($6 > '$maxNoCall' && NR != 1) print $1, $2 }' ${OUTPUT_DIR}common.hetnocall > ${OUTPUT_DIR}bad.common.txt
nBadRare=`wc -l ${OUTPUT_DIR}bad.rare.txt | awk '{print $1}'`
nBadCommon=`wc -l ${OUTPUT_DIR}bad.common.txt | awk '{print $1}'`
echo "Required the % of Het Calls for Rare SNPs to be between ${minHet}% and ${maxHet}% per sample ($nBadRare removed)" >> $LOG
cat ${OUTPUT_DIR}bad.rare.txt >> $LOG
echo "" >> $LOG
echo "Required the % of No Calls for Common SNPs to be less than ${maxNoCall}% per sample ($nBadCommon removed)" >> $LOG
cat ${OUTPUT_DIR}bad.common.txt >> $LOG

cat ${OUTPUT_DIR}bad.common.txt ${OUTPUT_DIR}bad.rare.txt > ${OUTPUT_DIR}bad.samples.txt

if [ ! -z $GAP_SUMMARY ]; then
    python ${scrd}addSummaryInfo.py ${OUTPUT_DIR}bad.samples.txt $GAP_SUMMARY > ${OUTPUT_DIR}badSampleInfo.txt
    echo "Counting the number of times a particular chip failed ..." >> $LOG
    awk '{print $3 }' ${OUTPUT_DIR}badSampleInfo.txt | sort -gk 1 | uniq -c >> $LOG
    python ${scrd}removeBadChips.py ${OUTPUT_DIR}badSampleInfo.txt $GAP_SUMMARY >> ${OUTPUT_DIR}bad.samples.txt
    echo "Removed samples for chips that had more than 6 samples failing sample QC..." >> $LOG   
fi

nBadTotal=`sort -gk 2 ${OUTPUT_DIR}bad.samples.txt | uniq -c | wc -l | awk '{print $1}'`
nRemain=$(($nSamples-$nBadTotal ))
echo "In total, removed $nBadTotal samples ... $nRemain samples remain" >> $LOG

echo "Samples Removed:" >> $LOG
cat ${OUTPUT_DIR}bad.samples.txt >> $LOG
        
date=`date +"%m/%d/%Y %T"`    
echo "Done with Sample QC -- $date" >> $LOG
for j in 1 2; do
    echo "" >> $LOG
done    

echo "$date" > ${OUTPUT_DIR}sampleqc.checksum
