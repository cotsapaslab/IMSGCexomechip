### USA-BSTN partition & MLMA workflow
### Mitrovic, Cotsapas Lab 2017


### Exercise 1: Partition the cohort based on PCA SNPs: 4 clusters (K = 4)
### Other options to test (--mcc, --cc, different # K)


## Define the command line argument
PLINK_file="$1"

## Partition into 4 clusters (toggle K for # of clusters)
./plink2 --bfile $PLINK_file --cluster --K 4 --out $PLINK_file
rm *.cluster3
rm *.cluster1

## Run PCA and produce first 20 PCs
./plink2 --bfile $PLINK_file --pca 20 header --out $PLINK_file

## Remove the header from the pca file
sed '1d' $PLINK_file.eigenvec > file1.temp

## Merge PCA (eigenvec) file with cluster solution file
paste file1.temp $PLINK_file.cluster2 | awk '{print $1,$2,$3,$4,$25}' > file2.temp

## Merge PC-cluster file with "SEX" and "PHENO" from fam file
paste file2.temp $PLINK_file.fam | awk '{print $1,$2,$3,$4,$5,$10,$11}' > final.table.temp

## Insert a header and produce a final file (7 columns)
awk 'BEGIN {print "FID SID PC1 PC2 CLUSTER PHENO SEX"} {print}' final.table.temp > final.table.K4.txt

## Check the number of cas/cons, males/females in each cluster
awk '{if ($5 == 0) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 1) print $2;}' | wc -l
awk '{if ($5 == 0) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 2) print $2;}' | wc -l
awk '{if ($5 == 1) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 1) print $2;}' | wc -l
awk '{if ($5 == 1) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 2) print $2;}' | wc -l
awk '{if ($5 == 2) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 1) print $2;}' | wc -l
awk '{if ($5 == 2) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 2) print $2;}' | wc -l
awk '{if ($5 == 3) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 1) print $2;}' | wc -l
awk '{if ($5 == 3) print $2,$6;}' final.table.K4.txt | awk '{if ($2 == 2) print $2;}' | wc -l

## Plot the clusters (saved in pca.plot.clusters.png)
Rscript pca_plot.clusters.R final.table.K4.txt

## Upper lines were compiled into the partition.sh script.


### Exercise 2: Partition the cohort based on PCA SNPs: 4 clusters and min 2000 cases and controls per cluster
### Other options to test (--cc, different # K)

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --mcc 2000 2000 --K 4 --out $PLINK_file


### Exercise 3: Partition the cohort based on PCA SNPs: 4 clusters and min 1000 cases and controls per cluster
### Other options to test (--cc, different # K)

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --mcc 1000 1000 --K 4 --out $PLINK_file


### Exercise 4: Partition the cohort based on PCA SNPs: 4 clusters and min 500 cases and controls per cluster
### Other options to test (--cc, different # K)

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --mcc 500 500 --K 4 --out $PLINK_file


### Exercise 5: Partition the cohort based on PCA SNPs: 3 clusters
### Other options to test (--cc, different # K)

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --K 3 --out $PLINK_file


### Exercise 6: Partition the cohort based on PCA SNPs: 3 clusters
### Other options to test (--cc, different # K)

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --K 3 --mcc 2000 2000 --out $PLINK_file

### Exercise 7: Partition the cohort based on PCA SNPs: 3 clusters
### Other options to test (--cc, different # K)

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --K 3 --mcc 1000 1000 --out $PLINK_file

### Exercise 8: Partition the cohort based on PCA SNPs: 5 clusters
### Other options to test (--cc, different # K)

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --K 5 --out $PLINK_file


### Exercise 9: Partition the cohort based on 121k MAF < 0.05 variants: 4 clusters
### Other options to test (--cc, different # K)

## PLINK =>  USA-BSTN.QCed.2nd_susp.noHLA.noXYMT.maf005

./plink2 --bfile $PLINK_file --cluster --K 4 --out $PLINK_file


### Exercise 10: Partition the cohort based on XXX MAF > 0.05 variants: 4 clusters
## Get the vars with MAF > 0.05
./plink2 --bfile USA-BSTN.QCed.2nd_susp --maf 0.05 --make-bed --out USA-BSTN.QCed.2nd_susp.minmaf005.temp1

## Exclude the MTXY vars
./plink2 --bfile USA-BSTN.QCed.2nd_susp.minmaf005.temp1 --not-chr 23, 24, 25, 26 --make-bed --out USA-BSTN.QCed.2nd_susp.minmaf005.temp2

## Exclude vars with missing call rate > 0.0001 (i. e. canonical call rate < 99.9999) 
./plink2 --bfile USA-BSTN.QCed.2nd_susp.minmaf005.temp2 --geno 0.0001 --make-bed --out USA-BSTN.QCed.2nd_susp.minmaf005.noMTXY

## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --K 4 --out $PLINK_file

### Exercise 11: Partition the cohort based on XXX MAF > 0.05 variants: 3 clusters
## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --K 3 --out $PLINK_file


### Exercise 12: Partition the cohort based on PCA SNPs: check the PCs 1,2,3,4
## Replace line 11 in partition.sh script with the following line and re-run the script
./plink2 --bfile $PLINK_file --cluster --K 4 --out $PLINK_file

sh partition.sh USA-BSTN.QCed.2nd_susp.pca.temp.pca.snps.temp

## Edited pca_plot.clusters.R to create pca_plot.clusters.v2.R and produce plots:
## - Plot PC1 vs PC2 by overlaying PHENO
## - Plot PC1 vs PC2 by overlaying SEX

## Edited partition.sh to create partition_v2.sh and produce plots:

## PC3 vs PC4 by PHENO
## PC1 vs PC3 by PHENO
## PC1 vs PC4 by PHENO
## PC2 vs PC3 by PHENO
## PC2 vs PC3 by PHENO
## PC2 vs PC4 by PHENO
## PC3 vs PC4 by SEX
## PC1 vs PC3 by SEX
## PC1 vs PC4 by SEX
## PC2 vs PC3 by SEX
## PC2 vs PC3 by SEX
## PC2 vs PC4 by SEX

## Edited pca_plot.clusters.v2.R to produce:
## Plot of BSTN cohorts as clusters on PC1 vs PC2 and PC3 vs PC4

### Exercise 13: Plot the 4 clusters on different PCs, remove outliers and run MLMA
## Remove cluster c02, c01 and c03 and check their sizes
awk '{if ($7 == 1) print $0}' final.table.K4.PCA_SNPs.txt > cluster.01.inds # 892 inds
awk '{if ($7 == 2) print $0}' final.table.K4.PCA_SNPs.txt > cluster.02.inds # 635 inds
awk '{if ($7 == 3) print $0}' final.table.K4.PCA_SNPs.txt > cluster.03.inds # 2 inds

## Remove inds from c1 - c03 from BSTN stratum
awk '{if ($7 != 0) print $0}' final.table.K4.PCA_SNPs.txt > clusters.01-03.inds # 1530 inds
./plink2 --bfile USA-BSTN.QCed.2nd_susp.pca.temp.pca.snps.temp --remove clusters.01-03.inds --make-bed --out USA-BSTN.QCed.2nd_susp.pca.temp.pca.snps.c0

## Re-run partiton_v2.sh on c0 inds of BSTN stratum (9466 inds, 5030 cases, 4435 controls, 2764 males, 6702 females)
sh partition_v2.sh USA-BSTN.QCed.2nd_susp.pca.temp.pca.snps.c0

## Plot clusters on PC plots using pca_plot.clusters.v3.R (manually, within R)

## Re-run MLMA on c0 inds of BSTN stratum (9466 inds, 5030 cases, 4435 controls, 2764 males, 6702 females)
## Extract c0 inds from the original BSTN stratum PLINK files
./plink2 --bfile USA-BSTN.QCed.2nd_susp --keep USA-BSTN.QCed.2nd_susp.pca.temp.pca.snps.c0.fam --make-bed --out USA-BSTN.QCed.2nd_susp.c0 
 
# Create phenotype covariate file (*.fam.phen) and sex covar file (*.fam.sex)
ls -1 USA-BSTN.QCed.2nd_susp.c0.fam > files.txt
for i in $(cat files.txt) ; do
awk '{print $1,$2,$6}' $i > $i'.phen'
awk '{print $1,$2,$5}' $i > $i'.sex'
done

###### NOTE: from here on, MLMA script runs on Farnam @ /ysm-gpfs/home/mm2534/BSTN.mlma/

# Create GRM
ls -1 USA-BSTN.QCed.2nd_susp.c0.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
./gcta64 --bfile $i --extract PCA_commonSNPs.autosomal.txt --make-bed --out $i'.pcaSNPs'
./gcta64 --bfile $i'.pcaSNPs' --make-grm-bin --out $i
done

# Run MLMA
ls -1 USA-BSTN.QCed.2nd_susp.c0.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
./gcta64 --mlma --bfile $i --grm-bin $i --pheno $i'.fam.phen' --covar $i'.fam.sex' --out $i
done

## Remove spaces from the MLMA assoc. file
ls -1 *.mlma > files.txt
for i in $(cat files.txt) ; do
sed 's/\ //g' $i > $i'.tab'
done

## Remove HLA variants from EMMAX assoc file
ls -1 *.mlma.tab > files.txt
for i in $(cat files.txt) ; do
awk -F " " 'BEGIN{while(getline<"hla.variants.extended.txt") a[$1]=1 } ; a[$2] !=1 {print $0 } ' $i > $i'.no.hla'
done

## Produce a Manhattan plot using manhattan.qq.R @ ~/Dropbox (Chris Cotsapas)/IMSGCechip/PCA/USA-BSTN.partition.and.mlma/

## Manually change header of USA-BSTN.QCed.2nd_susp.c0.results.mlma.tab.no.hla
## From Chr	SNP	bp	A1	A2	Freq	b	se	p
## To	CHR	SNP	BP	A1	A2	Freq	BETA	SE	P

## Run the script to produce the Manhattan and QQ plots

Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.c0.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.c0.fam


## Re-run MLMA on c1 and c2 inds of BSTN stratum in the same manner as for c0
./plink2 --bfile USA-BSTN.QCed.2nd_susp --keep cluster.01.inds --make-bed --out USA-BSTN.QCed.2nd_susp.c1 
./plink2 --bfile USA-BSTN.QCed.2nd_susp --keep cluster.02.inds --make-bed --out USA-BSTN.QCed.2nd_susp.c2


## Meta-analyze c0 + c1 + c2:
./plink2 --meta-analysis USA-BSTN.QCed.2nd_susp.c0.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.c1.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.c2.results.mlma.tab.no.hla + logscale --out USA-BSTN.QCed.2nd_susp.15kSNPs.c0.c1.c2

## Plot the meta results and get the lambdas:
Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.15kSNPs.c0.c1.c2.meta USA-BSTN.QCed.2nd_susp.fam

lambda SE lambda-1000
1.04966890657459 0.000185202090997923 1.00910246541845

### Exercise 14: Check if you get better clustering segregation with 22k home-made PCA SNPs
## Get the 22k SNPs
./plink2 --bfile USA-BSTN.QCed.2nd_susp --maf 0.05 --make-bed --out USA-BSTN.QCed.2nd_susp.minmaf005.temp1
## Exclude the MTXY vars
./plink2 --bfile USA-BSTN.QCed.2nd_susp.minmaf005.temp1 --not-chr 23, 24, 25, 26 --make-bed --out USA-BSTN.QCed.2nd_susp.minmaf005.temp2
## Exclude vars with missing call rate > 0.0001 (i. e. canonical call rate < 99.9999) 
./plink2 --bfile USA-BSTN.QCed.2nd_susp.minmaf005.temp2 --geno 0.0001 --make-bed --out USA-BSTN.QCed.2nd_susp.minmaf005.noMTXY

## Clustering and PCA
sh partition_v2.sh USA-BSTN.QCed.2nd_susp.minmaf005.noMTXY

awk '{if ($7 == 0) print $0}' final.table.K4.22k_SNPs.txt > cluster.22k_SNPs.c0.inds # 1919 inds (987 cases, 932 controls)
awk '{if ($7 == 1) print $0}' final.table.K4.22k_SNPs.txt > cluster.22k_SNPs.c1.inds # 1192 inds (695 cases, 497 controls)
awk '{if ($7 == 2) print $0}' final.table.K4.22k_SNPs.txt > cluster.22k_SNPs.c2.inds # 7883 inds (4286 cases, 3596 controls)
awk '{if ($7 == 3) print $0}' final.table.K4.22k_SNPs.txt > cluster.22k_SNPs.c3.inds # 1 con

./plink2 --bfile USA-BSTN.QCed.2nd_susp --keep cluster.22k_SNPs.c0.inds --make-bed --out USA-BSTN.QCed.2nd_susp.from22kSNPs.c0
./plink2 --bfile USA-BSTN.QCed.2nd_susp --keep cluster.22k_SNPs.c1.inds --make-bed --out USA-BSTN.QCed.2nd_susp.from22kSNPs.c1
./plink2 --bfile USA-BSTN.QCed.2nd_susp --keep cluster.22k_SNPs.c2.inds --make-bed --out USA-BSTN.QCed.2nd_susp.from22kSNPs.c2
./plink2 --bfile USA-BSTN.QCed.2nd_susp --keep cluster.22k_SNPs.c3.inds --make-bed --out USA-BSTN.QCed.2nd_susp.from22kSNPs.c3

###### NOTE: from here on, MLMA script runs on Farnam @ /ysm-gpfs/home/mm2534/BSTN.mlma/
sbatch run.mlma.22kSNPs.sh 

## Manually change header of *.results.mlma.tab.no.hla
## From Chr	SNP	bp	A1	A2	Freq	b	se	p
## To	CHR	SNP	BP	A1	A2	Freq	BETA	SE	P


## Run the script to produce the Manhattan and QQ plots

Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.from22kSNPs.c0.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.from22kSNPs.c0.fam
Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.from22kSNPs.c1.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.from22kSNPs.c1.fam
Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.from22kSNPs.c2.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.from22kSNPs.c2.fam

## Meta-analyze c0 + c1 + c2:
./plink2 --meta-analysis USA-BSTN.QCed.2nd_susp.from22kSNPs.c0.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.from22kSNPs.c1.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.from22kSNPs.c2.results.mlma.tab.no.hla + logscale --out USA-BSTN.QCed.2nd_susp.22kSNPs.c0.c1.c2

## Plot the meta results and get the lambdas:
Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.22kSNPs.c0.c1.c2.meta USA-BSTN.QCed.2nd_susp.fam



### Exercise 15: Check if you get no inflation without c2 or without c1 inds (derived from clustering on 15k SNPs)

## Remove the c2 (presumably Jewish pop) from the BSTN stratum
./plink2 --bfile USA-BSTN.QCed.2nd_susp --remove cluster.02.inds --make-bed --out USA-BSTN.QCed.2nd_susp.without.c2
./plink2 --bfile USA-BSTN.QCed.2nd_susp --remove cluster.01.inds --make-bed --out USA-BSTN.QCed.2nd_susp.without.c1

# Create phenotype covariate file (*.fam.phen) and sex covar file (*.fam.sex)
ls -1 USA-BSTN.QCed.2nd_susp.without.c1.fam > files.txt
for i in $(cat files.txt) ; do
awk '{print $1,$2,$6}' $i > $i'.phen'
awk '{print $1,$2,$5}' $i > $i'.sex'
done

ls -1 USA-BSTN.QCed.2nd_susp.without.c2.fam > files.txt
for i in $(cat files.txt) ; do
awk '{print $1,$2,$6}' $i > $i'.phen'
awk '{print $1,$2,$5}' $i > $i'.sex'
done

scp USA-BSTN.QCed.2nd_susp.without.c* mm2534@farnam.hpc.yale.edu:/ysm-gpfs/home/mm2534/BSTN.mlma/

###### NOTE: from here on, MLMA script runs on Farnam @ /ysm-gpfs/home/mm2534/BSTN.mlma/

sbatch run.mlma.without.c1.sh
sbatch run.mlma.without.c2.sh

###### NOTE: from here on we go back on local comp

scp mm2534@farnam.hpc.yale.edu:/ysm-gpfs/home/mm2534/BSTN.mlma/USA-BSTN.QCed.2nd_susp.without.c*.mlma* .

## Manually change header of *.results.mlma.tab.no.hla
## From Chr	SNP	bp	A1	A2	Freq	b	se	p
## To	CHR	SNP	BP	A1	A2	Freq	BETA	SE	P

## Plot the meta results and get the lambdas:
Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.without.c2.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.without.c2.fam
Rscript manhattan.qq.R USA-BSTN.QCed.2nd_susp.without.c1.results.mlma.tab.no.hla USA-BSTN.QCed.2nd_susp.without.c1.fam

## Get the first 20 PCs for the USA-BSTN.QCed.2nd_susp.without.c2 on 15k SNPs
## Run PCA and produce first 20 PCs
./plink2 --bfile USA-BSTN.QCed.2nd_susp.without.c2 --extract PCA_commonSNPs.autosomal.txt --make-bed --out USA-BSTN.QCed.2nd_susp.without.c2.PCA
./plink2 --bfile USA-BSTN.QCed.2nd_susp.without.c2.PCA --pca 20 header --out USA-BSTN.QCed.2nd_susp.without.c2.PCA

#### NOTE: The following command done here:
##/Users/mm/Dropbox (Chris Cotsapas)/IMSGCechip/PCA/USA-BSTN.partition.and.mlma/mlma.with.BSTN.partitioned

## MLMA meta with new BSTN stratum
./plink2 --meta-analysis Belgium.QCed.2nd_susp.mlma.tab.no.hla Greece.QCed.2nd_susp.mlma.tab.no.hla France.QCed.2nd_susp.mlma.tab.no.hla Finland.QCed.2nd_susp.mlma.tab.no.hla Netherlands.QCed.2nd_susp.mlma.tab.no.hla Norway.QCed.2nd_susp.mlma.tab.no.hla Italy.QCed.2nd_susp.mlma.tab.no.hla Denmark.QCed.2nd_susp.mlma.tab.no.hla Sweden.QCed.2nd_susp.mlma.tab.no.hla UK.Australia.QCed.2nd_susp.mlma.tab.no.hla USA-UCSF2.QCed.2nd_susp.mlma.tab.no.hla USA-UCSF1.QCed.2nd_susp.mlma.tab.no.hla Germany_eCHIP_LogReg_2PC_SEX.assoc.logistic.trim.tab.no.hla.betas USA-BSTN.QCed.2nd_susp.without.c2.results.mlma.tab.no.hla + logscale --out Final-metaanalysis-MLMA.no.hla.BSTN.wo.c2

## Create BSTN without c2 subgroup MLMA stats with p-vals and all counts 
## Working within /Users/mitja/Dropbox (Chris Cotsapas)/IMSGCechip/PCA/USA-BSTN.partition.and.mlma/ directory

## Get the counts (obtained by adjusted mac.sh)
./plink2 --bfile USA-BSTN.QCed.2nd_susp.without.c2 --freq --out USA-BSTN.QCed.2nd_susp.without.c2.maf --noweb
./plink2 --bfile USA-BSTN.QCed.2nd_susp.without.c2 --filter-cases --freq --counts --out USA-BSTN.QCed.2nd_susp.without.c2.counts.cases --noweb
./plink2 --bfile USA-BSTN.QCed.2nd_susp.without.c2 --filter-controls --freq --counts --out USA-BSTN.QCed.2nd_susp.without.c2.counts.controls --noweb
./plink2 --bfile USA-BSTN.QCed.2nd_susp.without.c2 --filter-cases --freq --out USA-BSTN.QCed.2nd_susp.without.c2.cases --noweb
./plink2 --bfile USA-BSTN.QCed.2nd_susp.without.c2 --filter-controls --freq --out USA-BSTN.QCed.2nd_susp.without.c2.controls --noweb
 
 paste USA-BSTN.QCed.2nd_susp.without.c2.maf.frq USA-BSTN.QCed.2nd_susp.without.c2.cases.frq USA-BSTN.QCed.2nd_susp.without.c2.controls.frq USA-BSTN.QCed.2nd_susp.without.c2.counts.cases.frq.counts USA-BSTN.QCed.2nd_susp.without.c2.counts.controls.frq.counts > USA-BSTN.QCed.2nd_susp.without.c2.temp
 awk '{print $1,$2,$3,$4,$5,$11,$17,$23,$24,$30,$31}' USA-BSTN.QCed.2nd_susp.without.c2.temp > USA-BSTN.QCed.2nd_susp.without.c2.final.counts

## Appropriate the header of the final counts file
 sed 's/.*CHR.*/CHR\ SNP\ A1\ A2\ MAF\ MAF_cas\ MAF_con\ Count_cas_A1 \Count_cas_A2 \ Count_con_A1 \ Count_con_A2/' USA-BSTN.QCed.2nd_susp.without.c2.final.counts > USA-BSTN.QCed.2nd_susp.without.c2.all.mlma.mafs.counts

## Match the MLMA stats with the counts by AWK oneliner
## MLMA (no HLA) has 147,038 lines
## Counts file has 150,887 lines

ls -1 USA-BSTN.QCed.2nd_susp.without.c2.results.mlma.tab.no.hla > files.txt
ls -1 USA-BSTN.QCed.2nd_susp.without.c2.all.mlma.mafs.counts > files1.txt
paste files.txt files1.txt |
while read i j; do
awk 'NR==FNR {h[$2] = $9; next} {print $0,h[$2]}' $i $j > $j'.with.pval'
done


## Always pass on what you have learned (Yoda)