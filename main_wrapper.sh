## This document describes detailed QC measures and analytical steps taken to produce the association
## results in the publication "Low frequency and rare coding variation contributes to multiple sclerosis risk".


## 1. Basic QC of cohorts
# 1.1 Run the QC pipeline:
bsub -o pipeline.out -R "rusage[mem=12]" "python /path/to/main.py path/to/meta"	 #LSF


## 2. Identity-by-missingness check (i.e. detection of possible correlations in patterns of missing data)

# 2.1 Calculate the IBM matrix
sh ibm.prep.sh

# 2.2 Calculate MDS components from IBM matrix, cluster on the resulting components, plot the results and output a list of outliers.
ls *QCed.autosomal.nohla.bad.out.fam > files.txt
for i in $(cat files.txt) ; do
 Rscript cluster.ibm.R $i
done

# 2.3 Remove potential IBM outliers
ls *.QCed.bed | sed -e 's/.bed$//g' > files1.txt
for i in $(cat files1.txt) ; do
	./plink --allow-no-sex \
	--bfile $i \
	--remove $i'.ibm.outliers.txt' \
	--make-bed \
	--out $i'.ibm.clean'
done


## 3. Projection PCA

# 3.1 Extract PCA SNPs
ls *.QCed.autosomal.nohla.bad.out.ibm.clean.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./plink --bfile $i \
	--extract PCA_commonSNPs.autosomal.txt \
	--make-bed \
	--out $i'.ibm.clean.pca.snps'
done

# 3.2 For optimal and fastest performance run pPCA script a cluster for each cohort
# Wrapper script and it's dependencies are in directory 3.pPCA
bsub -o pca.out -R "rusage[mem=8]" "sh PCA1_1KG_MSchip_exome_hg19.sh /path/to/1kg/ cohort1.QCed.autosomal.nohla.bad.out.ibm.clean.pca.snps EUR" #LSF

# 3.3 Remove possible pPCA outliers
ls *.QCed.ibm.clean.bed | sed -e 's/.bed$//g' > files1.txt
for i in $(cat files.txt) ; do
	./plink --allow-no-sex \
	--bfile $i \
	--remove $i'.pPCA.outliers.txt' \
	--make-bed \
	--out $i'.pca.clean'
done


## 4. Merging cohorts into strata

# 4.1 Before merging, updated allelic map to match European 1kG for each cohort.
# Main wrapper and it's dependencies are in directory 4.standardize_alleles.
sh standardize_alleles.echip.sh cohort1.QCed.ibm.clean.pca.clean.bim 

# 4.2 Cohorts were merged as shown in the strata.merging.scheme_cohorts.txt using PLINK, e.g.
./plink --bfile cohort1.QCed.ibm.clean.pca.clean --bmerge cohort2.QCed.ibm.clean.pca.clean.bed cohort2.QCed.ibm.clean.pca.clean.bim cohort2.QCed.ibm.clean.pca.clean.fam --make-bed --out stratum1


## 5. IBM on strata (merged cohorts)
# Repeat procedure from step 2 for each stratum.


## 6. Relatedness check in each stratum

# 6.1 Extract PCA SNPs for each stratum
ls *stratum.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
 	./plink --bfile $i \
	--extract PCA_commonSNPs.autosomal.txt \
 	--make-bed \
	--out $i'.ibd'\ 
 	--noweb
done

# 6.2 Calculate IBS matrix
ls *stratum.ibd.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./plink --bfile *.stratum.ibd \
	--genome \
	--out *.stratum.ibd
done

# 6.3 Calculate the missingness rates
ls *.stratum.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./plink --bfile  $i \ 
 	--missing \ 
 	--out $i \ 
 	--noweb
done

# 6.4 Identify possible IBD outliers (flag individual of the IBD pair that has the highest missing rate)
ls *stratum.ibd.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	perl run-IBD-QC.pl $i
done

# 6.5 Remove the IBD outliers
ls *stratum.bed | sed -e 's/.bed$//g' > files.txt
ls *stratum.ibd.fail-IBD-QC.txt > files1.txt
paste files.txt files1.txt |
while read i j; do
	./plink --bfile $i \ 
	--remove $j \ 
	--make-bed \
	--out $i'.ibdOK' \
	--noweb
done


## 7. Identify and remove variants with significant missingness differences between cases and controls

# 7.1 (Run the run-diff-miss.sh, which runs run-diffmiss-qc.pl in turn to get variants with significant missingness (p<0.0001) differences between cases and controls
ls *stratum.QCed.pca.clean.ibm.clean.ibdOK.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	sh run-diff-miss.sh 
done

## 8. pPCA for each stratum
# Repeat the procedure described in step 3.


## 9. Removal of monomorphic variants

# 9.1 Calculate MAF for all QCed variants for each stratum
ls *stratum.QCed.pca.clean.ibm.clean.ibdOK.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./plink --bfile $i \ 
	--freq \ 
	--out $i
done

# 9.2 Identify the monomorphic variants (i.e. those with MAF == 0)
ls *stratum.QCed.pca.clean.ibm.clean.ibdOK.frq
for i in $(cat files.txt) ; do
	awk '{if ($5 == "0") print $2;}' $i > $i'monomorphic.vars.txt'
done

# 9.3 Remove the identified monomorphic variants
ls *stratum.QCed.pca.clean.ibm.clean.ibdOK.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./plink --bfile $i \
	--exclude *$i'.monomorphic.vars.txt' \
	--make-bed \
	--out $i'.polymorphic'
done


# 10. Get the final PCs by flashPCA 
# flashPCA is a faster and equivalent alternative to EIGENSOFT

# 10.1 Calculate first 20 PCs
ls *stratum.QCed.pca.clean.ibm.clean.ibdOK.pca.snps.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./flashpca_x86-64 --bfile $i --ndim 20 --outpc $i'.20pcs'
done

# 10.2 Reformat and merge FID,SID,SEX,PHENO and PCs into one file
ls -1 *.pca.snps.fam > files1.txt
ls -1 *.20pcs > files2.txt
paste files1.txt files2.txt |
while read i j; do
	paste $i $j | cut -d' ' -f1,2,5,6,7- > $j'.temp'
done

ls -1 *.pca.snps.20pcs.temp > files.txt
for i in $(cat files.txt) ; do
	perl deeplink.pl $i tab
done


## 11. Mixed linear model association analysis (MLMA)

# 11.1 Create phenotype covariate file (*.fam.phen) and sex covar file (*.fam.sex)
ls -1 *.QCed.stratum.fam > files.txt
for i in $(cat files.txt) ; do
	awk '{print $1,$2,$6}' $i > $i'.phen'
	awk '{print $1,$2,$5}' $i > $i'.sex'
done

# 11.2 Create GRM
ls -1 *.QCed.stratum.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./gcta64 --bfile $i --extract PCA_commonSNPs.autosomal.txt --make-bed --out $i'.pcaSNPs'
	./gcta64 --bfile $i'.pcaSNPs' --make-grm-bin --out $i
done

# 11.3 Run MLMA
ls -1 *.QCed.stratum.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	bsub -o mlma.out -R "rusage[mem=16]" "./gcta64 --mlma --bfile $i --grm-bin $i --pheno $i'.fam.phen' --covar $i'.fam.sex' --out $i" #LSF
done

## 11.4 Remove the spaces from the MLMA assoc. file
ls -1 *.mlma > files.txt
	for i in $(cat files.txt) ; do
	sed 's/\ //g' $i > $i'.tab'
done

## 11.5 Remove HLA variants from EMMAX assoc file
ls -1 *.mlma.tab > files.txt
for i in $(cat files.txt) ; do
	awk -F " " 'BEGIN{while(getline<"hla.variants.extended.txt") a[$1]=1 } ; a[$2] !=1 {print $0 } ' $i > $i'.no.hla'
done

## 12. Meta-analysis of per stratum association statistics

./plink --meta-analysis Belgium.mlma.tab.no.hla Denmark.mlma.tab.no.hla Greece.mlma.tab.no.hla France.mlma.tab.no.hla Finland.mlma.tab.no.hla Sweden.mlma.tab.no.hla Norway.mlma.tab.no.hla Netherlands.mlma.tab.no.hla Italy.mlma.tab.no.hla UK.Australia.mlma.tab.no.hla USA-BSTN.mlma.tab.no.hla USA-UCSF1.mlma.tab.no.hla USA-UCSF2.mlma.tab.no.hla Germany_eCHIP_LogReg_2PC_SEX.assoc.logistic.trim.tab.no.hla.betas + logscale --out echip.mlma.meta.no.hla




