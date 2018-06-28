# IMSGC ExomeChip QC and analysis pipeline
QC and analysis code used to produce the association results published in "Low frequency and rare coding variation contributes to multiple sclerosis risk"

IMSGC exome chip QC/analysis pipeline <br>
International Multiple Sclerosis Genetics Consortium <br>
21 June 2018 <br>
[Mitja Mitrovic](mailto:mitja.mitrovic@yale.edu) <br>
[Chris Cotsapas](mailto:cotsapas@broadinstitute.org) <br>

## Front matter
This README file documents the quality control (QC) and analysis code used to produce the association results published in [Low frequency and rare coding variation contributes to multiple sclerosis risk](https://www.biorxiv.org/content/early/2018/03/23/286617). Here we provide an overview of the software provided in this distribution, the software dependencies, and the sources of genotype data required. We also provide lists of samples and SNPs that pass our QC process, if you do not want to step through the entire QC pipeline.  

## Disclaimers
* Throughout, we use the LSF job scheduling system to distribute some of the more compute-intensive analyses on a cluster. We have tagged all instances with `#LSF` in `Main_wrapper.txt` for ease of reference.
* Association statistics calculated with different numbers of samples and/or stratum structure will vary from those we report. 
* The code here is provided as-is. We have validated that this runs without errors on our own systems, but do not guarantee that this will be the case in different computing environments.

## Genotype data acquisition and structure
Pre-QC genotype data for 36,219 cases and 38,629 controls can be obtained from the IMSGC (dac@imsgc.net). 

Our samples are split into 42 cohorts which represent sample collections, systematic batches and chip versions. In the paper, we assembled our cohorts into 13 strata for analysis. For each cohort, we provide:
+ A README file
+ Raw genotypes from Gencall, as plink bed/bim/fam files
+ Genotypes after zCall processing of Gencall files, as plink bed/bim/fam files
+ A list of samples passing our QC, in plink fam format, named `<prefix>.qcpass.inds`
+ A list of variants passing our QC, in plink bim format, named `<prefix>.qcpass.snps`
+ A data processing masterfile as an xlsx Excel workbook, showing how many samples and variants pass each step of our QC pipeline for all cohorts in the collection.
+ A list of all cohorts belonging to the same stratum as the current cohort, as a text file of plink file name roots. The file is named `<prefix>-merge-scheme.txt`
+ A lookup table for variant IDs used in these files to canonical rsIDs (hg19/gb137). NA means no rsID record exists for that variant in hg19/gb137. The file is named `20180402-rsIDs.txt`


## Recreating our post-QC dataset 
You can recreate our dataset without running raw data through this pipeline. To do this, extract the individuals and markers passing QC in each cohort (per-cohort `*.exc.qcpass.inds and *.exc.qcpass.snps` files). You can then choose to recreate our analysis strata (per-cohort `<prefix>-merge-scheme.txt` files).

## Running our QC pipeline
Here, we provide the full code required to re-analyze our raw data and replicate our published results.  

### Pipeline overview
Our pipeline is split into twelve distinct steps, including meta-analysis of association results across strata (summarized in Figure S1 of our manuscript). These are:

+ Basic QC per cohort: remove markers and samples with high missingness, Hardy-Weinberg filters, heterozygosity by missingness filter. 
+ Remove samples affected by batch effects per cohort (identity-by-missingness filter)
+ Remove population outliers per cohort by projection PCA
+ Merge cohorts into strata
+ Remove samples affected by batch effects per stratum (identity-by-missingness filter)
+ Remove related individuals per stratum
+ Remove markers with differential missingness between cases and controls
+ Remove population outliers per stratum by projection PCA
+ Remove monomorphic variants
+ Calculate final PCA for population stratification control 
+ Calculate case/control association statistics with mixed linear models per stratum
+ Meta-analyze association statistics across all strata

### Software dependencies
+ GNU bash 3.2.57 or later
+ EIGENSOFT 6.1.3 or later
+ GCTA v1.91.4
+ FlashPCA 1.2.6 or later
+ Perl v5.12 or later
+ PLINK v1.9
+ Python 3.0
+ R v2.15 or later (including packages ggplot2 v2.2.1, mclust v5.4, GenABEL v1.8, qqman v0.1.4, dplyr v0.7.4)

### Running the pipeline
The entire pipeline can be run from `main_wrapper.sh`, which calls the scripts included in this distribution to execute all twelve steps of our pipeline. There are detailed instructions for each step in the wrapper script. User input is required at various stages:

+ In Step 1, for each cohort you must do …. 
+ In Step 1, for each thing you must do …
+ In Step 5, you must have the HLA variant file in the same directory … 

### Notes about options
In step 12, note that the `+ logscale` option is required because the per-stratum association files list effect betas, NOT odds ratios (remember that `beta = ln(OR)` ).


### Structure and description of files and directories in this distribution

+ Directory/file name
+ Description
1 basic_qc_pipeline
QC scripts and suppl. files.
2 post_qc_ibm
Identity-by-missingness scripts and suppl. files.
3 pPCA
Projection PCA scripts and suppl. files.
4 standardize_alleles
Scripts and suppl. files to update allelic map to match European 1000 Genomes.
7 differential_missingness
Two scripts to identify SNPs with differentially missing data between cases and controls.
Main_wrapper.sh
Wrapper with 12 subsequent steps from raw files gencall and zcall files to MLMA meta-analysis.
PCA_commonSNPs.autosomal.txt
A list of common, LD-pruned autosomal SNPs used for PCA. 
Per.cohort.minHet.maxHet.xlsx
Min. and max. heterozygosity cutoffs for each cohort, used in meta file of QC pipeline.
README.txt
This file.
strata.merging.scheme_cohorts.txt
A file showing how the 42 cohorts were merged into 13 strata.


1.basic_qc_pipeline directory

This directory contains scripts and other supplementary files needed to run the basic QC pipeline. Detailed instructions, pipeline workflow and options are described in ~/code4github/1.basic_qc_pipeline/README. Before running the pipeline, as described in step 1.1 in the Main_wrapper.sh, you have to modify ~/code4github/1.basic_qc_pipeline/meta file. Example structure of the meta file used in this analysis is described in ~/code4github/1.basic_qc_pipeline/README and shown in ~/code4github/1.basic_qc_pipeline/meta.

Per cohort minimal and maximal heterozygosity used in this QC are saved in ~/code4github/Per.cohort.minHet.maxHet.xlsx. When processing a dataset with N > 3,000, it is recommended to run the QC on cluster with a week wall-time.


2.post_qc_ibm

Steps 2.1-3 in the Main_wrapper.sh provide the commands for running the IBM check. Script ibm.prep.sh first extracts non-HLA autosomal variants, which are then used in  calculation of identity-by-missingness matrix. Please note that ibm.prep.sh scripts expects both lists of variants, sex.mt.chr.variants.txt and hla.variants.extended.txt, to be present in the same directory. Next, cluster.ibm.R is used to calculate the MDS components from IBM matrix to detect possible correlations in patterns of missing data and output a list of outliers. Finally, potentially identified outliers are removed step 2.3.


3.pPCA

Steps 3.1-3 in the Main_wrapper.sh show the commands for running pPCA. The code expects the list of common, autosomal and LD-pruned variants saved under PCA_commonSNPs.autosomal.txt in the ~/code4github/3.pPCA/ directory. It is suggested that the script that executes pPCA, PCA1_1KG_MSchip_exome_hg19.sh, is run on a cluster for optimal performance. Please note that the second argument in the line 45 of the Main_wrapper.sh, /path/to/1kg/, is a path to PLINK level European 1000 Genomes data (phase 3) parsed per chromosome (e.g. chr1.bed/bim/fam, chr2.bed/bim/fam etc). These files are not provided in this distribution. Projection PCA is performed on the input data and on input data merged with 1000 Genomes individuals. Outlier detection is based on 6 standard deviations (smartpca default setting). Outliers are removed within 10 iterations (smartpca default setting), where each iteration consists of following steps: calculation of PCs, detection of outliers (individuals outside 6 SDs), removal of outliers, re-calculation of PCs.


4.standardize_alleles

Scripts in this directory update allelic maps of all cohorts to match the European (EUR) 1000 Genomes  before merging them into strata. The main script, standardize_alleles.echip.sh, expects the following dependencies in the same directory: 
*  check_alleles.sh, 
*  check_alleles.header, 
*  allele_map_comp_flipped_to_original.dat, 
*  allele_map_comp_to_original.dat, 
*  allele_map_flipped_to_original.dat, 
*  echip_original_alleles.dat.

Once on the same allelic map, cohorts are then merged into the 13 strata according to scheme described in the ~/code4github/strata.merging.scheme_cohorts.txt and command in step 4.2 in the Main_wrapper.sh.

6. Relatedness check in each stratum

The code in step 6.1 of the Main_wrapper.sh expects the list of common, autosomal and LD-pruned variants saved under PCA_commonSNPs.autosomal.txt in the same directory.

Step 6.4 of the Main_wrapper.sh expects the script run-IBD-QC.pl in the same directory.


7.differential_missingness

The main script expects run-diffmiss-qc.pl in the same directory. Significance threshold is set to p < 0.0001 by default and can be changed in line 14 of the run-diffmiss-qc.pl script.

11. MLMA

Step 11.2 of the Main_wrapper.sh expects PCA_commonSNPs.autosomal.txt in the same directory. Since MLMA is computationally intensive, it is suggested to run it on cluster, as shown for LSF in line 199.

Step 11.5 of the Main_wrapper.sh expects a list of HLA variants, hla.variants.extended.txt, in the same directory.


12. Meta-analysis of 13 strata MLMA association statistics

Please note that “+ logscale” option was used since the variant effects are already on log-scale (i.e. betas from MLMA).
