#!/bin/bash

## Paths to pipeline directories
src_direc=/IMSGCexomechip/src
bin_direc=/IMSGCexomechip/bin
temp_direc=/IMSGCexomechip/temp
strata_direc=/IMSGCexomechip/temp/strata
suppl_direc=/IMSGCexomechip/suppl


PATH=$PATH:${bin_direc}

################################################################################
############################    Step 0: Notes       ############################
################################################################################
## This script runs raw ExomeChip PLINK files through entire pipeline including 
## detailed QC measures and analytical steps to produce the association results 
## in our publication "Low frequency and rare coding variation contributes to 
## multiple sclerosis risk".
##
## Our pipeline is split into twelve distinct steps:
##
## 1.  Basic QC per cohort: remove markers and samples with high missingness, 
## 	   Hardy-Weinberg filters, heterozygosity by missingness filter. 
## 2.  Remove samples affected by batch effects per cohort
## 3.  Remove population outliers per cohort by projection PCA
## 4.  Merge cohorts into strata
## 5.  Remove samples affected by batch effects per stratum
## 6.  Remove related individuals per stratum
## 7.  Remove markers with differential missingness between cases and controls
## 8.  Remove population outliers per stratum by projection PCA
## 9.  Remove monomorphic variants
## 10. Calculate final PCA for population stratification control 
## 11. Calculate case/control association statistics with mixed linear models 
##     per stratum
## 12. Meta-analyze association statistics across all strata
##
## This script is a wrapper, which calls the scripts included in this distribution
## to execute all twelve steps of our pipeline. User input is required at various 
## stages as described in detail in the README.txt.


################################################################################
############################    Step 1: Basic QC    ############################
################################################################################

# 1.1 Modify ${src_direc}/basic_qc/meta file

# 1.2 Run the QC pipeline for ExomeChip cohorts on a computer cluster:
#LSF	
bsub -o pipeline.out -R "rusage[mem=12]" "python ${src_direc}/basic_qc/main.py ${src_direc}/basic_qc/meta"

################################################################################
###################   Step 2: Identity-by-missingness check   ##################
################################################################################

# 2.1 List all QCed cohorts
ls -1 ${temp_direc}/basic_qc_out/ > cohorts.txt

# 2.2 Calculate the IBM matrix per cohort
for cohort in $(cat cohorts.txt); do
	mkdir -p ${temp_direc}/ibm/${cohort}/
	mkdir -p ${temp_direc}/logs/${cohort}/
	
	./plink --bfile ${temp_direc}/basic_qc_out/${cohort}/${cohort} \
			--exclude ${suppl_direc}/sex.mt.chr.variants.txt \
			--make-bed \
			--allow-no-sex \
			--out ${temp_direc}/ibm/${cohort}/${cohort}.autosomal

	./plink --bfile ${temp_direc}/ibm/${cohort}.autosomal \
			--exclude ${suppl_direc}/hla.variants.extended.txt \
			--make-bed \
			--allow-no-sex \
			--out ${temp_direc}/ibm/${cohort}/${cohort}.autosomal.nohla
	#LSF	
	bsub -o ${temp_direc}/logs/${cohort}/ibm.prep.log.out.${cohort}. "./plink --bfile ${temp_direc}/ibm/${cohort}/${cohort}.autosomal.nohla --cluster-missing --allow-no-sex --out ${temp_direc}/ibm/${cohort}/${cohort}.autosomal.nohla"
done	
	
# 2.2 Calculate MDS components from IBM matrix, cluster on the resulting 
# components, plot the results and output a list of outliers.
for cohort in $(cat cohorts.txt); do
 Rscript ${src_direc}/cluster.ibm.R ${temp_direc}/ibm/${cohort}/${cohort}.autosomal.nohla
done

# 2.3 Remove IBM outliers
for cohort in $(cat cohorts.txt); do
	./plink --bfile ${temp_direc}/basic_qc_out/${cohort}/${cohort} \
			--remove ${temp_direc}/ibm/${cohort}/${cohort}.autosomal.nohla.ibm.outliers.txt \
			--make-bed \
			--allow-no-sex \
			--out ${temp_direc}/ibm/${cohort}/${cohort}.ibm.clean
done

################################################################################
###########################  Step 3: Projection PCA  ###########################
################################################################################

# 3.1 Extract PCA SNPs
for cohort in $(cat cohorts.txt); do
	mkdir -p ${temp_direc}/pPCA/${cohort}/
	./plink --bfile ${temp_direc}/ibm/${cohort}/${cohort}.ibm.clean \
			--extract ${suppl_direc}/PCA_commonSNPs.autosomal.txt \
			--make-bed \
			--allow-no-sex \
			--out ${temp_direc}/pPCA/${cohort}/${cohort}.pca.snps
done

# 3.2 Run pPCA on a computer cluster for each cohort
#LSF	
bsub -o ${temp_direc}/logs/${cohort}/ppca.log.out.${cohort} "sh ${src_direc}/pPCA/PCA1_1KG_MSchip_exome_hg19.sh /path/to/1kg/ ${temp_direc}/pPCA/${cohort}/${cohort}.pca.snps EUR"

# 3.3 Remove pPCA outliers
for cohort in $(cat cohorts.txt); do
	./plink --allow-no-sex \
			--bfile ${temp_direc}/ibm/${cohort}/${cohort}.ibm.clean \
			--remove ${temp_direc}/pPCA/${cohort}/${cohort}.pPCA.outliers.txt \
			--make-bed \
			--out ${temp_direc}/pPCA/${cohort}/${cohort}.pca.clean
done


################################################################################
########################  Step 4: Merging cohorts into strata  #################
################################################################################

# 4.1 Before merging, updated allelic map to match European 1kG for each cohort
for cohort in $(cat cohorts.txt); do
	mkdir -p ${temp_direc}/merge/${cohort}/
	cd ${temp_direc}/merge/${cohort}/
	sh ${src_direc}/standardize_alleles.echip.sh ${temp_direc}/pPCA/${cohort}/${cohort}.pca.clean
done

# 4.2 Cohorts were merged as shown in ${suppl_direc}/strata.merging.scheme_cohorts.txt
# Example below shows merging of Danish cohorts into "Denmark" stratum:
./plink --bfile DEN1.ibm.clean.pca.clean.stand_alleles \
		--merge-list Danish.cohorts.merge \
		--make-bed
		--allow-no-sex \
		--out Denmark
		
#cat Danish.cohorts.merge
#DenmarkPt1.ibm.clean.pca.clean.stand_alleles.bed DenmarkPt1.ibm.clean.pca.clean.stand_alleles.bim cDenmarkPt1.ibm.clean.pca.clean.stand_alleles.fam
#DenmarkPt2.ibm.clean.pca.clean.stand_alleles.bed DenmarkPt2.ibm.clean.pca.clean.stand_alleles.bim cDenmarkPt2.ibm.clean.pca.clean.stand_alleles.fam

# 4.3 At this point we copy binary PLINK level strata files to ${strata_direc}


################################################################################
###################   Step 5: Identity-by-missingness check   ##################
################################################################################


# 5.1 Calculate the IBM matrix per cohort
for stratum in $(cat ${suppl_direc}/strata.txt); do
	mkdir -p ${strata_direc}/ibm/${stratum}
	
	./plink --bfile ${strata_direc} \
			--exclude ${suppl_direc}/sex.mt.chr.variants.txt \
			--make-bed \
			--allow-no-sex \
			--out ${strata_direc}/ibm/${stratum}.autosomal

	./plink --bfile ${strata_direc}/ibm/${stratum}/${stratum}.autosomal \
			--exclude ${suppl_direc}/hla.variants.extended.txt \
			--make-bed \
			--allow-no-sex \
			--out ${strata_direc}/ibm/${stratum}/${stratum}.autosomal.nohla
	
	#LSF	
	bsub -o ${strata_direc}/ibm/${stratum}/ibm.log.out.${stratum} "./plink --bfile ${strata_direc}/ibm/${stratum}/${stratum}.autosomal.nohla --cluster-missing --allow-no-sex --out ${strata_direc}/ibm/${stratum}/${stratum}.autosomal.nohla"
done	
	
# 5.2 Calculate MDS components from IBM matrix, cluster on the resulting 
# components, plot the results and output a list of outliers.
for stratum in $(cat ${suppl_direc}/strata.txt); do
 Rscript ${src_direc}/cluster.ibm.R ${strata_direc}/ibm/${stratum}/${stratum}.autosomal.nohla
done

# 5.3 Remove IBM outliers
for stratum in $(cat ${suppl_direc}/strata.txt); do
	./plink --bfile ${strata_direc}/ibm/${stratum}/${stratum} \
			--remove ${strata_direc}/ibm/${stratum}/${stratum}.autosomal.nohla.ibm.outliers.txt \
			--make-bed \
			--allow-no-sex \
			--out ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean
done


################################################################################
########################   Step 6: Relatedness check  ##########################
################################################################################


# 6.1 Extract PCA SNPs for each stratum and calculate IBS matrix
for stratum in $(cat ${suppl_direc}/strata.txt); do
	mkdir -p ${strata_direc}/ibd/${stratum}
 	./plink --bfile ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean \
			--extract ${suppl_direc}/PCA_commonSNPs.autosomal.txt \
 			--make-bed \
			--out ${strata_direc}/ibd/${stratum}/${stratum}.autosomes \ 
	#LSF	
	bsub -o ${strata_direc}/ibd/${stratum}/ibd.log.out.${stratum} "./plink --bfile ${strata_direc}/ibd/${stratum}/${stratum}.autosomes --genome --allow-no-sex --out ${strata_direc}/ibd/${stratum}/${stratum}.autosomes"
done

# 6.2 Calculate the missingness rates per stratum, then identify possible IBD outliers 
# (flag individual of the IBD pair that has the highest missing rate) and finally,
# remove the IBD outliers
for stratum in $(cat ${suppl_direc}/strata.txt); do
	./plink --bfile  ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean \
 			--missing \ 
 			--out ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean \ 
 			--allow-no-sex

	perl ${src_direc}/run-IBD-QC.pl ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean

	./plink --bfile ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean \ 
			--remove ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean.fail-IBD-QC.txt \ 
			--make-bed \
			--out ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean.ibd.clean \
			--allow-no-sex
done


################################################################################
###################   Step 7: Differential missingness check  #################
################################################################################


# 7.1 Identify and remove variants with significant missingness differences between cases and controls
for stratum in $(cat ${suppl_direc}/strata.txt); do	
	mkdir -p ${strata_direc}/diffmiss/${stratum}/
	./plink --bfile ${strata_direc}/ibm/${stratum}/${stratum}.ibm.clean.ibd.clean \
			--test-missing \
			--out ${strata_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean \
			--allow-no-sex

	perl ${src_direc}/run-diffmiss-qc.pl ${strata_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean
	#echo ${strata_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean.fail-diffmiss-qc.txt | \
	#awk '{split($0,a,"."); print a[1]}'
	#cat $i | wc -l
	./plink --bfile ${strata_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean \
			--exclude ${temp_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean.fail-diffmiss-qc.txt \
			--make-bed \
			--out ${strata_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean \
			--allow-no-sex
done


################################################################################
######################   Step 8: Projection PCA per stratum ####################
################################################################################

# 8.1 Extract PCA SNPs
for stratum in $(cat ${suppl_direc}/strata.txt); do
	mkdir -p ${strata_direc}/ppca/${stratum}/
	./plink --bfile ${strata_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean \
			--extract ${suppl_direc}/PCA_commonSNPs.autosomal.txt \
			--make-bed \
			--allow-no-sex \
			--out ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.snps
done

# 8.2 Run pPCA on a computer cluster for each stratum
#LSF	
bsub -o ${strata_direc}/ppca/${stratum}/ppca.log.out.${stratum} "sh ${src_direc}/pPCA/PCA1_1KG_MSchip_exome_hg19.sh /path/to/1kg/ ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.snps EUR"

# 8.3 Remove pPCA outliers
for stratum in $(cat ${suppl_direc}/strata.txt); do
	./plink --allow-no-sex \
			--bfile ${strata_direc}/diffmiss/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean \
			--remove ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.snps.pPCA.outliers.txt \
			--make-bed \
			--out ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.clean
done

################################################################################
###################   Step 9: Remove monomorphic SNPs per stratum ##############
################################################################################

# 9.1 Calculate MAF, identify the monomorphic variants (i.e. those with MAF == 0) and
# remove the identified monomorphic variants
mkdir -p ${strata_direc}/QCed.strata/
for stratum in $(cat ${suppl_direc}/strata.txt); do
	./plink --bfile ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.clean \ 
			--freq \
			--allow-no-sex \
			--out ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.clean
	awk '{if ($5 == "0") print $2;}' ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.clean.frq > \
	${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.clean.monomorphic.vars.txt

	./plink --bfile ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.clean \
	--exclude ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.clean.monomorphic.vars.txt \
	--make-bed \
	--allow-no-sex \
	--out ${strata_direc}/QCed.strata/${stratum}.QCed
done

################################################################################
###################   Step 10: Obtain final PCs per stratum ####################
################################################################################

# 10.1 Calculate first 20 PCs
for stratum in $(cat ${suppl_direc}/strata.txt); do
	./flashpca_x86-64 --bfile ${strata_direc}/ppca/${stratum}/${stratum}.ibm.clean.ibd.clean.diff-miss.clean.pca.snps --ndim 20 --outpc ${strata_direc}/QCed.strata/${stratum}.QCed.20pcs
done


################################################################################
##############   Step 11: Mixed linear model association analysis (MLMA) #######
################################################################################

# 11.1 Create phenotype covariate file (*.fam.phen) and sex covar file (*.fam.sex)
# and create GRM
for stratum in $(cat ${suppl_direc}/strata.txt); do
	mkdir -p ${strata_direc}/mlma/${stratum}/
	awk '{print $1,$2,$6}' ${strata_direc}/QCed.strata/${stratum}.QCed.fam > ${strata_direc}/QCed.strata/${stratum}.QCed.phen
	awk '{print $1,$2,$5}' ${strata_direc}/QCed.strata/${stratum}.QCed.fam > ${strata_direc}/QCed.strata/${stratum}.QCed.sex
	
	./gcta64 --bfile ${strata_direc}/QCed.strata/${stratum}.QCed \
			 --extract ${suppl_direc}/PCA_commonSNPs.autosomal.txt \
			 --make-bed \
			 --out ${strata_direc}/mlma/${stratum}.QCed.pcaSNPs
	#LSF
	bsub -o ${strata_direc}/mlma/${stratum}.grm.log.out.${stratum} "./gcta64 --bfile ${strata_direc}/mlma/${stratum}.QCed.pcaSNPs --make-grm-bin --out ${strata_direc}/mlma/${stratum}.QCed"
done

# 11.2 Run MLMA and produce the assoc. statistics
	#LSF
for stratum in $(cat ${suppl_direc}/strata.txt); do
	bsub -o o ${strata_direc}/mlma/${stratum}.mlma.log.out.${stratum} -R "./gcta64 --mlma --bfile --bfile ${strata_direc}/mlma/${stratum}.QCed --grm-bin ${strata_direc}/mlma/${stratum}.QCed --pheno ${strata_direc}/QCed.strata/${stratum}.QCed.phen --covar ${strata_direc}/QCed.strata/${stratum}.QCed.sex --out ${strata_direc}/QCed.strata/${stratum}.QCed"
done

## 11.3 Remove space characters and HLA variants from MLMA assoc. file
for stratum in $(cat ${suppl_direc}/strata.txt); do
	sed 's/\ //g' ${strata_direc}/mlma/${stratum}.QCed.mlma > ${strata_direc}/mlma/${stratum}.QCed.mlma.tab
	awk -F " " 'BEGIN{while(getline<"${suppl_direc}/hla.variants.extended.txt") a[$1]=1 } ; a[$2] !=1 {print $0 } ' ${strata_direc}/mlma/${stratum}.QCed.mlma.tab > ${strata_direc}/mlma/${stratum}.QCed.mlma.tab.no.HLA
done

################################################################################
##############  Step 12:Meta-analyze association statistics across all strata ##
################################################################################

# 12.1 MLMA meta-analysis with HLA region variants
./plink --meta-analysis Belgium.QCed.mlma.tab Denmark.QCed.mlma.tab Greece.QCed.mlma.tab France.QCed.mlma.tab Finland.QCed.mlma.tab Sweden.QCed.mlma.tab Norway.QCed.mlma.tab Netherlands.QCed.mlma.tab Italy.QCed.mlma.tab UK.Australia.QCed.mlma.tab USA-BSTN.QCed.mlma.tab USA-UCSF1.QCed.mlma.tab USA-UCSF2.QCed.mlma.tab Germany_eCHIP_LogReg_2PC_SEX.assoc.logistic.trim.tab.no.hla.betas + logscale --out echip.mlma.meta

# 12.2 MLMA meta-analysis without HLA region variants
./plink --meta-analysis Belgium.QCed.mlma.tab.no.hla Denmark.QCed.mlma.tab.no.hla Greece.QCed.mlma.tab.no.hla France.QCed.mlma.tab.no.hla Finland.QCed.mlma.tab.no.hla Sweden.QCed.mlma.tab.no.hla Norway.QCed.mlma.tab.no.hla Netherlands.QCed.mlma.tab.no.hla Italy.QCed.mlma.tab.no.hla UK.Australia.QCed.mlma.tab.no.hla USA-BSTN.QCed.mlma.tab.no.hla USA-UCSF1.QCed.mlma.tab.no.hla USA-UCSF2.QCed.mlma.tab.no.hla Germany_eCHIP_LogReg_2PC_SEX.assoc.logistic.trim.tab.no.hla.betas + logscale --out echip.mlma.meta.no.hla


