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

## Define variables with directories
wd=`grep "Script Directory:" $META | awk '{print $3}'`
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/


date=`date +"%m/%d/%Y %T"`
echo "[$i] PCA using Eigenstrat -- $date" >> $LOG; ((i++))

pcSites=${supd}pca_snps.txt        
plink --bfile $GENCALL --remove ${OUTPUT_DIR}bad.samples.txt --geno $GC_CR --make-bed --out ${OUTPUT_DIR}gencall.tmp1 --extract $pcSites --maf 0.05

python ${scrd}findDupsWControls.py ${OUTPUT_DIR}gencall.tmp1.fam ${cond}controls.fam > ${OUTPUT_DIR}overlapWControls.txt
q=`wc -l ${OUTPUT_DIR}overlapWControls.txt | awk '{ print $1 }'`
if [ $(($q)) -ne 0 ]; then
    plink --bfile ${OUTPUT_DIR}gencall.tmp1 --remove ${OUTPUT_DIR}overlapWControls.txt --make-bed --out ${OUTPUT_DIR}tmp
    mv ${OUTPUT_DIR}tmp.bed ${OUTPUT_DIR}gencall.tmp1.bed
    mv ${OUTPUT_DIR}tmp.bim ${OUTPUT_DIR}gencall.tmp1.bim
    mv ${OUTPUT_DIR}tmp.fam ${OUTPUT_DIR}gencall.tmp1.fam                        
fi
        
COHORT=`grep "Cohort File:" $META | awk '{print $3}'`
if [ -z "$COHORT" ]; then            
    echo "FID IID cohort" > ${OUTPUT_DIR}cohort.txt
    awk -v pn="$PROJECT_NAME" '{ print $1, $2, pn }' $GENCALL.fam >> ${OUTPUT_DIR}cohort.txt
    COHORT=${OUTPUT_DIR}cohort.txt
fi

## Convertf Data
python ${scrd}convertPhenoToPedInd.py ${OUTPUT_DIR}gencall.tmp1.fam $COHORT > ${OUTPUT_DIR}gencall.tmp1.pedind

echo "genotypename:    ${OUTPUT_DIR}gencall.tmp1.bed" > ${OUTPUT_DIR}par.eigenstrat
echo "snpname:         ${OUTPUT_DIR}gencall.tmp1.bim # or example.map, either works" >> ${OUTPUT_DIR}par.eigenstrat
echo "indivname:       ${OUTPUT_DIR}gencall.tmp1.pedind # or example.ped, either works" >> ${OUTPUT_DIR}par.eigenstrat
echo "outputformat:    EIGENSTRAT" >> ${OUTPUT_DIR}par.eigenstrat
echo "genotypeoutname: ${OUTPUT_DIR}tmp.geno" >> ${OUTPUT_DIR}par.eigenstrat
echo "snpoutname:      ${OUTPUT_DIR}tmp.snp" >> ${OUTPUT_DIR}par.eigenstrat
echo "indivoutname:    ${OUTPUT_DIR}tmp.ind" >> ${OUTPUT_DIR}par.eigenstrat
echo "familynames:     YES    " >> ${OUTPUT_DIR}par.eigenstrat

${eigd}bin/convertf -p ${OUTPUT_DIR}par.eigenstrat > ${OUTPUT_DIR}convertf.log
echo "Finished converting data to eigenstrat format..." >> $LOG

## Merge with 1kg/controls data
echo "geno1: ${cond}controlsPCA.geno" > ${OUTPUT_DIR}par.mergeit
echo "snp1:  ${cond}controlsPCA.snp" >> ${OUTPUT_DIR}par.mergeit
echo "ind1:  ${cond}controlsPCA.ind" >> ${OUTPUT_DIR}par.mergeit
echo "geno2: ${OUTPUT_DIR}tmp.geno" >> ${OUTPUT_DIR}par.mergeit
echo "snp2:  ${OUTPUT_DIR}tmp.snp" >> ${OUTPUT_DIR}par.mergeit
echo "ind2:  ${OUTPUT_DIR}tmp.ind" >> ${OUTPUT_DIR}par.mergeit
echo "genotypeoutname: ${OUTPUT_DIR}merge.geno" >> ${OUTPUT_DIR}par.mergeit
echo "snpoutname:      ${OUTPUT_DIR}merge.snp" >> ${OUTPUT_DIR}par.mergeit
echo "indivoutname:    ${OUTPUT_DIR}merge.ind" >> ${OUTPUT_DIR}par.mergeit
echo "outputformat: EIGENSTRAT" >> ${OUTPUT_DIR}par.mergeit
${eigd}bin/mergeit -p ${OUTPUT_DIR}par.mergeit > ${OUTPUT_DIR}mergeit.log
echo "Merged with controls data..." >> $LOG
    
## smartPCA with outlier detection
echo "genotypename:    ${OUTPUT_DIR}merge.geno" > ${OUTPUT_DIR}par.smartpca
echo "snpname:         ${OUTPUT_DIR}merge.snp" >> ${OUTPUT_DIR}par.smartpca
echo "indivname:       ${OUTPUT_DIR}merge.ind" >> ${OUTPUT_DIR}par.smartpca
echo "evecoutname:     ${OUTPUT_DIR}merge.evec" >> ${OUTPUT_DIR}par.smartpca
echo "evaloutname:     ${OUTPUT_DIR}merge.eval" >> ${OUTPUT_DIR}par.smartpca
echo "altnormstyle:    NO # use price normalization" >> ${OUTPUT_DIR}par.smartpca
echo "numoutevec:      10" >> ${OUTPUT_DIR}par.smartpca
echo "numoutlieriter: 0" >> ${OUTPUT_DIR}par.smartpca
echo "poplistname: ${supd}populationListForEigenstrat.txt" >> ${OUTPUT_DIR}par.smartpca
    
${eigd}bin/smartpca -p ${OUTPUT_DIR}par.smartpca > ${OUTPUT_DIR}smartpca.log
Rscript ${scrd}plotEVEC.r ${OUTPUT_DIR}merge.evec ${OUTPUT_DIR}pcaround1.evec.plot.pc12.png ${OUTPUT_DIR}pcaround1.evec.plot.first4PCs.png ${OUTPUT_DIR}pcaround1.evec.ind.plot.png
        
## Convert smartPCA eigenvectors into covariates file
python ${scrd}makeCovar.py ${OUTPUT_DIR}merge.evec > ${OUTPUT_DIR}covar.all.txt
## Identify European subset        
awk '{ if ($3 > 0.006) print $1, $2 }' ${OUTPUT_DIR}covar.all.txt > ${OUTPUT_DIR}european.txt
## Identify African subset        
awk '{ if ($4 > 0.02) print $1, $2 }' ${OUTPUT_DIR}covar.all.txt > ${OUTPUT_DIR}african.txt
## Identify Asian subset        
awk '{ if ($4 < -0.027 && $3 < -0.02) print $1, $2 }' ${OUTPUT_DIR}covar.all.txt > ${OUTPUT_DIR}asian.txt

nEur=`python ${scrd}countEthnSamples.py ${GENCALL}.fam ${OUTPUT_DIR}european.txt | awk '{print $1}'`
nAfr=`python ${scrd}countEthnSamples.py ${GENCALL}.fam ${OUTPUT_DIR}african.txt | awk '{print $1}'`
nAsian=`python ${scrd}countEthnSamples.py ${GENCALL}.fam ${OUTPUT_DIR}asian.txt | awk '{print $1}'`        
        
echo "Finished subsetting samples into different populations..." >> $LOG    
echo "European: PC1 > 0.006 ... $nEur samples" >> $LOG
echo "African: PC2 > 0.02 ... $nAfr samples" >> $LOG
echo "Asian: PC1 < -0.02 and PC2 < -0.027 ... $nAsian samples" >> $LOG
echo "Ignoring samples from the Americas..." >> $LOG
    
for j in 1 2; do
    echo "" >> $LOG
done

echo "$date" > ${OUTPUT_DIR}pcaround1.checksum