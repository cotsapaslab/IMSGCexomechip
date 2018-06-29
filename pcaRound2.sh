#!/bin/bash

alias plink='plink --noweb --silent --allow-no-sex'

## Define variables
OUTPUT_DIR=$2
COHORT=$4
PROJECT_NAME=$3

## Define variables with directories
wd=$1
scrd=${wd}/scripts/
eigd=${wd}/eigenstrat/
hlad=${wd}/hlaImputation/
cond=${wd}/controls/
supd=${wd}/suppl/
        
pcSites=${supd}pca_snps.txt        
plink --bfile ${OUTPUT_DIR}${PROJECT_NAME}_CONTROLS --make-bed --out ${OUTPUT_DIR}pckeep --extract $pcSites --keep ${OUTPUT_DIR}european.txt --remove ${cond}nonEuropean.txt
                
## Convertf Data
python ${scrd}convertPhenoToPedInd.py ${OUTPUT_DIR}pckeep.fam $COHORT > ${OUTPUT_DIR}pckeep.pedind

echo "genotypename:    ${OUTPUT_DIR}pckeep.bed" > ${OUTPUT_DIR}par.eigenstrat
echo "snpname:         ${OUTPUT_DIR}pckeep.bim" >> ${OUTPUT_DIR}par.eigenstrat
echo "indivname:       ${OUTPUT_DIR}pckeep.pedind" >> ${OUTPUT_DIR}par.eigenstrat
echo "outputformat:    EIGENSTRAT" >> ${OUTPUT_DIR}par.eigenstrat
echo "genotypeoutname: ${OUTPUT_DIR}pckeep.geno" >> ${OUTPUT_DIR}par.eigenstrat
echo "snpoutname:      ${OUTPUT_DIR}pckeep.snp" >> ${OUTPUT_DIR}par.eigenstrat
echo "indivoutname:    ${OUTPUT_DIR}pckeep.ind" >> ${OUTPUT_DIR}par.eigenstrat
echo "familynames:     YES    " >> ${OUTPUT_DIR}par.eigenstrat

${eigd}bin/convertf -p ${OUTPUT_DIR}par.eigenstrat > ${OUTPUT_DIR}convertf.log


## smartPCA with outlier detection
echo "genotypename:    ${OUTPUT_DIR}pckeep.geno" > ${OUTPUT_DIR}par.smartpca
echo "snpname:         ${OUTPUT_DIR}pckeep.snp" >> ${OUTPUT_DIR}par.smartpca
echo "indivname:       ${OUTPUT_DIR}pckeep.ind" >> ${OUTPUT_DIR}par.smartpca
echo "evecoutname:     ${OUTPUT_DIR}euro.evec" >> ${OUTPUT_DIR}par.smartpca
echo "evaloutname:     ${OUTPUT_DIR}euro.eval" >> ${OUTPUT_DIR}par.smartpca
echo "altnormstyle:    NO # use price normalization" >> ${OUTPUT_DIR}par.smartpca
echo "numoutevec:      10" >> ${OUTPUT_DIR}par.smartpca
echo "numoutlieriter: 0" >> ${OUTPUT_DIR}par.smartpca
echo "poplistname: ${supd}populationListForEigenstrat.txt" >> ${OUTPUT_DIR}par.smartpca
    
${eigd}bin/smartpca -p ${OUTPUT_DIR}par.smartpca > ${OUTPUT_DIR}smartpca.log
Rscript ${scrd}plotEVEC2.r ${OUTPUT_DIR}euro.evec ${OUTPUT_DIR}euro.evec.plot.pc12.png ${OUTPUT_DIR}euro.evec.plot.first4PCs.png ${OUTPUT_DIR}euro.evec.ind.plot.png

## Convert smartPCA eigenvectors into covariates file
python ${scrd}makeCovar.py ${OUTPUT_DIR}euro.evec > ${OUTPUT_DIR}covar.euro.txt

echo "" > ${OUTPUT_DIR}pcaround2.checksum
