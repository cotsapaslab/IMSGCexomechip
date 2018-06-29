#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Define variables
OUTPUT_DIR=$2
PROJECT_NAME=$3

## Define variables with directories
wd=$1
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

maxJobs=$4

pcSites=${supd}pca_snps.txt
plink --bfile ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS --make-bed --out ${OUTPUT_DIR}pcKeep0 --extract $pcSites --maf 0.05
plink --bfile ${OUTPUT_DIR}pcKeep0 --freq --out ${OUTPUT_DIR}tmp.freq

if [ ! -d ${OUTPUT_DIR}bsubLogs ]; then
    mkdir ${OUTPUT_DIR}bsubLogs
else
    a=`ls -l ${OUTPUT_DIR}bsubLogs/ | wc -l | awk '{print $1 }'`
    if [ $(($a)) -ne 0 ]; then
        rm ${OUTPUT_DIR}bsubLogs/*
    fi
fi

gawk '{print $1,$2}' ${OUTPUT_DIR}pcKeep0.fam | split -d -a 3 -l 500 - ${OUTPUT_DIR}tmp.list

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
        bsub -P ${PROJECT_NAME} -q hour -o ${OUTPUT_DIR}bsubLogs/$z.log plink --bfile ${OUTPUT_DIR}pcKeep0 --read-freq ${OUTPUT_DIR}tmp.freq.frq --genome --genome-lists ${OUTPUT_DIR}tmp.list`printf "%03i\n" $i` ${OUTPUT_DIR}tmp.list`printf "%03i\n" $j` --out ${OUTPUT_DIR}data.sub.$i.$j --noweb --allow-no-sex
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

cat ${OUTPUT_DIR}data.sub*genome | awk '{ if (NR == 1 || (!(match($0, "FID")))) print $0 }' > ${OUTPUT_DIR}genome.genome
rm ${OUTPUT_DIR}tmp.list*
rm ${OUTPUT_DIR}data.sub.*

awk '{ if (NR == 1 || $10 > 0.2) print $0 }' ${OUTPUT_DIR}genome.genome > ${OUTPUT_DIR}genome_0.2.genome                        
awk '{ if ($10 > 0.7 && NR != 1) print $0 }' ${OUTPUT_DIR}genome_0.2.genome > ${OUTPUT_DIR}dups.txt
awk '{ if ($8 > 0.7 && NR != 1) print $0 }' ${OUTPUT_DIR}genome_0.2.genome > ${OUTPUT_DIR}offspring.txt
awk '{ if (($7 > 0.2 && $8 > 0.4 && $9 > 0.2) && NR != 1) print $0 }' ${OUTPUT_DIR}genome_0.2.genome > ${OUTPUT_DIR}sibs.txt            

nDups=`wc -l ${OUTPUT_DIR}dups.txt | awk '{print $1}'`

if [ $(($nDups)) -gt 0 ]; then
    awk '{print $1, $2}' ${OUTPUT_DIR}dups.txt > ${OUTPUT_DIR}tmp.dups.txt
    awk '{print $3, $4}' ${OUTPUT_DIR}dups.txt >> ${OUTPUT_DIR}tmp.dups.txt                
    awk '{print $3, $4 }' ${OUTPUT_DIR}dups.txt > ${OUTPUT_DIR}bad.samples.txt
    cat ${OUTPUT_DIR}tmp.dups.txt | sort -gk 1 | uniq -c | awk '{ if ($1 > 1) print $2, $3 }' >> ${OUTPUT_DIR}bad.samples.txt    
    rm ${OUTPUT_DIR}tmp.dups.txt

    plink --bfile ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS --remove ${OUTPUT_DIR}bad.samples.txt --make-bed --out ${OUTPUT_DIR}tmp

    mv ${OUTPUT_DIR}tmp.fam ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.fam
    mv ${OUTPUT_DIR}tmp.bed ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.bed
    mv ${OUTPUT_DIR}tmp.bim ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.bim
    mv ${OUTPUT_DIR}tmp.log ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS.log    
fi
                                                                             
echo "$date" > ${OUTPUT_DIR}relatedness_check2.checksum
