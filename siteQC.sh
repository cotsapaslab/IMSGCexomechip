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

GC_CR=`grep "gencall_missrate:" $META | awk '{print $2}'`
ZC_CR=`grep "zcall_missrate:" $META | awk '{print $2}'`

## Define variables with directories
wd=`grep "Script Directory:" $META | awk '{print $3}'`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/

knownBadSNPs=${supd}badSNPs.txt

## Check for Errors
function die(){
    echo "Error in siteQC.sh..."
    exit 1
    }
#trap die ERR

date=`date +"%m/%d/%Y %T"`
echo "[$i] Site QC -- $date" >> $LOG; ((i++))
echo "Using GenCall Maximum No Call Rate of $GC_CR" >> $LOG
echo "Using zCall Maximum No Call Rate of $ZC_CR" >> $LOG        

python ${scrd}removeSibsOffspring.py ${OUTPUT_DIR}sibs.txt ${OUTPUT_DIR}offspring.txt > ${OUTPUT_DIR}tmp.bad.samples.txt
cat ${OUTPUT_DIR}bad.samples.txt >> ${OUTPUT_DIR}tmp.bad.samples.txt #List of unrelated samples

if [ -e ${OUTPUT_DIR}sexCheckErrors.txt ]; then
    cat ${OUTPUT_DIR}sexCheckErrors.txt >> ${OUTPUT_DIR}tmp.bad.samples.txt
fi
if [ -e ${OUTPUT_DIR}mendelCheckErrors.txt ]; then
    cat ${OUTPUT_DIR}mendelCheckErrors.txt >> ${OUTPUT_DIR}tmp.bad.samples.txt
fi

awk '{ if (NR != 1) print $1,$2 }' ${OUTPUT_DIR}european.txt > ${OUTPUT_DIR}keepEuro.txt

## QC Autosomal SNPs        
plink --bfile $GENCALL --remove ${OUTPUT_DIR}tmp.bad.samples.txt --geno $GC_CR --make-bed --out ${OUTPUT_DIR}gencall.tmp2
awk '{ if ($1 < 23 || $1 == 25) print $2 }' ${OUTPUT_DIR}gencall.tmp2.bim  > ${OUTPUT_DIR}keepSNPS.gcAuto.txt
plink --bfile $ZCALL --remove ${OUTPUT_DIR}tmp.bad.samples.txt --keep ${OUTPUT_DIR}keepEuro.txt --geno $ZC_CR --hwe 1e-5 --make-bed --out ${OUTPUT_DIR}zCall.autosomal.eur --extract ${OUTPUT_DIR}keepSNPS.gcAuto.txt
awk '{print $2 }' ${OUTPUT_DIR}zCall.autosomal.eur.bim > ${OUTPUT_DIR}keepAutosomalSNPs.txt
plink --bfile $ZCALL --remove ${OUTPUT_DIR}bad.samples.txt --geno $ZC_CR --make-bed --out ${OUTPUT_DIR}zCall.qc.auto --extract ${OUTPUT_DIR}keepAutosomalSNPs.txt --exclude $knownBadSNPs
nAuto=`wc -l ${OUTPUT_DIR}zCall.qc.auto.bim | awk '{print $1 }'`
echo "Finished QC'ing Autosomal and Pseudo-Autosomal SNPs ... kept $nAuto snps" >> $LOG    

## QC X Chr SNPs (initial)
awk '{ if ($1 == 23) print $2 }' ${OUTPUT_DIR}gencall.tmp2.bim  > ${OUTPUT_DIR}keepSNPS.gcX.txt    
plink --bfile $ZCALL --remove ${OUTPUT_DIR}tmp.bad.samples.txt --keep ${OUTPUT_DIR}keepEuro.txt --filter-females --hwe 1e-5 --geno $ZC_CR --make-bed --out ${OUTPUT_DIR}zcall.tmp3 --extract ${OUTPUT_DIR}keepSNPS.gcX.txt --update-sex ${OUTPUT_DIR}sexUpdate.txt
awk '{ print $2 }' ${OUTPUT_DIR}zcall.tmp3.bim > ${OUTPUT_DIR}keepX.txt
plink --bfile $ZCALL --remove ${OUTPUT_DIR}bad.samples.txt --extract ${OUTPUT_DIR}keepX.txt --geno $ZC_CR --make-bed --out ${OUTPUT_DIR}zCall.qc.x    
plink --bfile ${OUTPUT_DIR}zCall.qc.auto --bmerge ${OUTPUT_DIR}zCall.qc.x.bed ${OUTPUT_DIR}zCall.qc.x.bim ${OUTPUT_DIR}zCall.qc.x.fam --make-bed --out ${OUTPUT_DIR}zCall.qc
nX=`wc -l ${OUTPUT_DIR}zCall.qc.x.bim | awk '{print $1 }'`
echo "Finished QC'ing X SNPs ... kept $nX snps" >> $LOG    

## QC Y Chr SNPs
awk '{ if ($1 == 24) print $2 }' $GENCALL.bim > ${OUTPUT_DIR}keepY.txt
plink --bfile $GENCALL --remove ${OUTPUT_DIR}tmp.bad.samples.txt --extract ${OUTPUT_DIR}keepY.txt --filter-males --geno $GC_CR --recode --out ${OUTPUT_DIR}y --transpose --update-sex ${OUTPUT_DIR}sexUpdate.txt
python ${scrd}countHetTPED.py ${OUTPUT_DIR}y.tped > ${OUTPUT_DIR}y.hetnocall
awk '{ if ($2 < 3) print $1 }' ${OUTPUT_DIR}y.hetnocall > ${OUTPUT_DIR}keepY.txt
plink --bfile $ZCALL --remove ${OUTPUT_DIR}bad.samples.txt --extract ${OUTPUT_DIR}keepY.txt --make-bed --out ${OUTPUT_DIR}zCall.qc.y
nY=`wc -l ${OUTPUT_DIR}zCall.qc.y.bim | awk '{print $1 }'`
echo "Finished QC'ing Y SNPs ... kept $nY snps" >> $LOG    

## Merge Good Autosomal, X, Y SNPs back together
plink --bfile ${OUTPUT_DIR}zCall.qc --bmerge ${OUTPUT_DIR}zCall.qc.y.bed ${OUTPUT_DIR}zCall.qc.y.bim ${OUTPUT_DIR}zCall.qc.y.fam --make-bed --out ${OUTPUT_DIR}${PROJECT_NAME}

echo "Finished QC'ing SNPs" >> $LOG    
for j in 1 2; do
    echo "" >> $LOG
done

echo "$date" > ${OUTPUT_DIR}snpqc.checksum