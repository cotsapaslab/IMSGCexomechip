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
maxJobs=`grep "maxJobs:" $META | awk '{ print $2 }'`

## Define variables with directories
wd=`grep "Script Directory:" $META | awk '{print $3}'`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/


GC_CR=`grep "gencall_missrate:" $META | awk '{print $2}'`

## Check for Errors
function die(){
    echo "Error in relatednessCheck.sh..." >> $LOG
    exit 1
    }
#trap die ERR

date=`date +"%m/%d/%Y %T"`
echo "[$i] Begin Relatedness Check -- $date" >> $LOG; ((i++))

pcSites=${supd}pca_snps.txt
plink --bfile $GENCALL --remove ${OUTPUT_DIR}bad.samples.txt --geno $GC_CR --make-bed --out ${OUTPUT_DIR}gencall.pcKeep0 --extract $pcSites --maf 0.05
plink --bfile ${OUTPUT_DIR}gencall.pcKeep0 --freq --out ${OUTPUT_DIR}tmp.freq
if [ ! -d ${OUTPUT_DIR}bsubLogs ]; then
    mkdir ${OUTPUT_DIR}bsubLogs
else
    a=`ls -l ${OUTPUT_DIR}bsubLogs/ | wc -l | awk '{print $1 }'`
    if [ $(($a)) -ne 0 ]; then
        rm ${OUTPUT_DIR}bsubLogs/*
    fi
fi

echo "Splitting up genome calculation in PLINK" >> $LOG
gawk '{print $1,$2}' ${OUTPUT_DIR}gencall.pcKeep0.fam | split -d -a 3 -l 500 - ${OUTPUT_DIR}tmp.list

i=0
a=`ls -l ${OUTPUT_DIR}tmp.list* | wc -l | awk '{print $1 }'`
a=$(($a - 1))
j=0
z=0
while [ $i -le $a ]; do
    while [ $j -le $a ]; do
        nJobs=`bjobs | wc -l`
        if [ $(($nJobs)) -gt $(($maxJobs)) ]; then
            sleep 5m
            continue
        fi    
        ((z++))
        bsub -P ${PROJECT_NAME} -q hour -o ${OUTPUT_DIR}bsubLogs/$z.log plink --bfile ${OUTPUT_DIR}gencall.pcKeep0 --read-freq ${OUTPUT_DIR}tmp.freq.frq --genome --genome-lists ${OUTPUT_DIR}tmp.list`printf "%03i\n" $i` ${OUTPUT_DIR}tmp.list`printf "%03i\n" $j` --out ${OUTPUT_DIR}data.sub.$i.$j --noweb --allow-no-sex
        ((j++))
    done
    ((i++))
    j=$i
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
cat ${OUTPUT_DIR}data.sub*genome | awk '{ if (NR == 1 || (!(match($0, "FID")))) print $0 }' > ${OUTPUT_DIR}genome.genome
rm ${OUTPUT_DIR}tmp.list*
rm ${OUTPUT_DIR}data.sub.*
            
echo "Looking for dups, parent/offspring, siblings, and overly related samples..." >> $LOG

awk '{ if (NR == 1 || $10 > 0.2) print $0 }' ${OUTPUT_DIR}genome.genome > ${OUTPUT_DIR}genome_0.2.genome
awk '{ if (NR == 1) print $1,$2,$3,$4,$7,$8,$9,$10 }' ${OUTPUT_DIR}genome_0.2.genome >> $LOG
awk '{ if ($10 > 0.7 && NR != 1) print $0 }' ${OUTPUT_DIR}genome_0.2.genome > ${OUTPUT_DIR}dups.txt
nDups=`wc -l ${OUTPUT_DIR}dups.txt | awk '{print $1}'`
if [ $(($nDups)) -ne 0 ]; then
    awk '{print $1, $2, $3, $4, $7,$8,$9, $10}' ${OUTPUT_DIR}dups.txt >> $LOG
fi
            
awk '{ if ($8 > 0.7 && NR != 1) print $0 }' ${OUTPUT_DIR}genome_0.2.genome > ${OUTPUT_DIR}offspring.txt
nOffspring=`wc -l ${OUTPUT_DIR}offspring.txt | awk '{print $1}'`
if [ $(($nOffspring)) -ne 0 ]; then
    awk '{print $1, $2, $3, $4, $7,$8,$9, $10}' ${OUTPUT_DIR}offspring.txt >> $LOG
fi

awk '{ if (($7 > 0.2 && $8 > 0.4 && $9 > 0.2) && NR != 1) print $0 }' ${OUTPUT_DIR}genome_0.2.genome > ${OUTPUT_DIR}sibs.txt            
nSibs=`wc -l ${OUTPUT_DIR}sibs.txt | awk '{print $1}'`
if [ $(($nSibs)) -ne 0 ]; then
    awk '{print $1, $2, $3, $4, $7,$8,$9, $10}' ${OUTPUT_DIR}sibs.txt >> $LOG
fi
            
python ${scrd}countPiHatAboveThresh.py ${OUTPUT_DIR}genome.genome | awk '{ if (NR > 1 && $3 > 5) print $0 }' > ${OUTPUT_DIR}piHatCount.txt  
nPiHat=`wc -l ${OUTPUT_DIR}piHatCount.txt | awk '{print $1}'`

echo "" >> $LOG

STUDY_TYPE=`grep "Study Type:" $META | awk '{print $3}'`                        
exit=0
if [ $(($nDups)) -gt 0 ]; then
    echo "Duplicates exist in file. Automatically removing FID2,IID2. If this is not the desired result, remove desired samples from input files and restart pipeline." >> $LOG
    awk '{print $1,$2, $3, $4 }' ${OUTPUT_DIR}dups.txt >> $LOG
    awk '{print $1, $2}' ${OUTPUT_DIR}dups.txt > ${OUTPUT_DIR}tmp.dups.txt
    awk '{print $3, $4}' ${OUTPUT_DIR}dups.txt >> ${OUTPUT_DIR}tmp.dups.txt                
    awk '{print $3, $4 }' ${OUTPUT_DIR}dups.txt >> ${OUTPUT_DIR}bad.samples.txt
    cat ${OUTPUT_DIR}tmp.dups.txt | sort -gk 1 | uniq -c | awk '{ if ($1 > 1) print $2, $3 }' >> ${OUTPUT_DIR}bad.samples.txt
    rm ${OUTPUT_DIR}tmp.dups.txt
    echo "" >> $LOG                
fi

## if parent/offspring undefined in plink files -- > warning to redefine in fam file --> exit 1
if [ $(($nOffspring)) -gt 0 ]; then
    python ${scrd}checkParOffspring.py ${OUTPUT_DIR}offspring.txt $GENCALL.fam > ${OUTPUT_DIR}parentOffspringCheck.txt
    nFalse=`wc -l ${OUTPUT_DIR}parentOffspringCheck.txt | awk '{print $1}'`
    if [ $nFalse -gt 0 ]; then
        echo "Parent/Offspring relationships not defined in fam file." >> $LOG
        awk '{print $1, $2, $3, $4 }' ${OUTPUT_DIR}parentOffspringCheck.txt >> $LOG
        remove=`grep "unrelated_output:" $META | awk '{ print $2 }'`
        if [ $remove = 1 ]; then
            awk '{print $3, $4 }' ${OUTPUT_DIR}parentOffspringCheck.txt >> ${OUTPUT_DIR}bad.samples.txt
        fi
        echo "" >> $LOG                
    fi
fi

if [ $(($nSibs)) -gt 0 ]; then
    python ${scrd}checkSibs.py ${OUTPUT_DIR}sibs.${n}.txt $GENCALL.fam > ${OUTPUT_DIR}sibCheck.txt
    nFalse=`wc -l ${OUTPUT_DIR}sibCheck.txt | awk '{print $1}'`
    if [ $nFalse -gt 0 ]; then
        echo "Sibling relationships not defined in fam file." >> $LOG
        awk '{print $1, $2, $3, $4 }' ${OUTPUT_DIR}sibCheck.txt >> $LOG
        remove=`grep "unrelated_output:" $META | awk '{ print $2 }'`
        if [ $remove = 1 ]; then            
            awk '{print $3, $4 }' ${OUTPUT_DIR}sibCheck.txt >> ${OUTPUT_DIR}bad.samples.txt
        fi                    
        echo "" >> $LOG                    
    fi
fi

if [ $(($nPiHat)) -gt 0 ]; then
    echo "Identified $nPiHat samples that are overly related to other samples. If samples are all from one population (ex: european), then consider removing samples in [ ${OUTPUT_DIR}piHatCount.$n.txt ] and restarting the pipeline. " >> $LOG
    echo "" >> $LOG                
fi


date=`date +"%m/%d/%Y %T"`    
echo "Done with Relatedness Check -- $date" >> $LOG

if [ $exit -eq 1 ]; then
    exit 1
fi

for j in 1 2; do
    echo "" >> $LOG
done    

echo "$date" > ${OUTPUT_DIR}relatedness_check.checksum