# IMSGC ExomeChip QC and analysis pipeline

## Contact
International Multiple Sclerosis Genetics Consortium <br>
21 June 2018 <br>
[Mitja Mitrovic](mailto:mitja.mitrovic@yale.edu) <br>
[Chris Cotsapas](mailto:cotsapas@broadinstitute.org) <br>

## Front matter
This README file documents the quality control (QC) and analysis code used to produce the association results published in [Low frequency and rare coding variation contributes to multiple sclerosis risk](https://www.biorxiv.org/content/early/2018/03/23/286617). Here we provide an overview of the software provided in this distribution, the software dependencies, and the sources of genotype data required. We also provide lists of samples and SNPs that pass our QC process, if you do not want to step through the entire QC pipeline.  

## Disclaimers
* Throughout, we use the LSF job scheduling system to distribute some of the more compute-intensive analyses on a cluster. We have tagged all instances with `#LSF` in `main.pipeline.sh` for ease of reference.
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
The entire pipeline can be run from `main.pipeline.sh`, which calls the scripts included in this distribution to execute all twelve steps of our pipeline. There are detailed instructions for each step in the script. The script is not intended as a stand-alone, unsupervised piece of code - rather, it steps through the various procedures in our QC approach. User input is required at various stages:

+ In step 1, you must first modify the meta file. An example meta file used in our analysis is in IMSGCexomechip/src/basic_qc/meta. Structure of meta file used   in this analysis is described in ~/IMSGCexomechip/1.basic_qc_pipeline/README. You must provide path to your input zcall and gencall PLINK files and output directory. In addition, copy controls_bu.zip and controls.zip from our DropBox [link] and expand them to the /IMSGCexomechip/src/basic_qc/ before running the pipeline.
+ In step 2.2, careful inspection of cluster plots before outlier removal is advised.
+ In step 3.2, before running the PCA1_1KG_MSchip_exome_hg19.sh script on computer cluster, you must specify path to your copy of 1000 Genomes data (Phase 3), which should be in binary (bed/bim/fam) PLINK format and parsed per chromosome (e.g. chr1.bed/bim/fam, chr2.bed/bim/fam etc). These files are not provided in this distribution.
+ In step 8.2, before running the PCA1_1KG_MSchip_exome_hg19.sh script on computer cluster, you must specify path to your copy of 1000 Genomes data (Phase 3), which should be in binary (bed/bim/fam) PLINK format and parsed per chromosome (e.g. chr1.bed/bim/fam, chr2.bed/bim/fam etc). These files are not provided in this distribution.
+ In step 10, we use FlashPCA, as a faster and equivalent alternative to EIGENSOFT.


### Notes about options
+ In step 1, detailed instructions, pipeline workflow and options are described in /IMSGCexomechip/src/basic_qc/README.txt 
+ In step 1.2, when processing datasets with N > 3,000 individuals, it is recommended to run the QC on a computer cluster with a week wall-time.
+ In step 3.2, the third variable is the population under study, which was EUR in our experiment. All possible populations are: EUR ASN AFR AMR. Outlier detection is based on 6 standard deviations (smartpca default setting). Outliers are removed within 10 iterations (smartpca default setting), where each iteration consists of following steps: calculation of PCs, detection of outliers (individuals outside 6 SDs), removal of outliers, re-calculation of PCs.
+ In step 6.1, when processing datasets with N > 3,000 individuals, cluster analyses can be very slow. We recommend using a compute cluster to run these calculations in parallel, as described in PLINK notes.
+ In step 7.1, which runs run-diffmiss-qc.pl to identify variants with significant missingness differences between cases and controls, the default p is set to p<0.0001.
+ In step 10, we use FlashPCA to calculate the first twenty principal components by default. This can be adjusted with the --ndim option.
+ In step 12, note that the `+ logscale` option is required because the per-stratum association files list effect betas, NOT odds ratios (remember that `beta = ln(OR)` ).


### Description of files and directories in this distribution

| File/directory | Description |
| ------------- | ------------- |
| bin/  | Contains executable versions of PLINK, GCTA and FlashPCA.   |
| src/  | Contains all scripts required in by the pipeline.  |
| main.pipeline.sh  | Main script with twelve QC and analysis steps of our pipeline.  |
| README.txt  | This file.  |
| supp/hla.variants.extended.txt  | List of 4,629 ExomeChip HLA region variants, chr6:24,046,865-34,985,625. |
| supp/PCA_commonSNPs.autosomal.txt  | List of 16,066 common (MAF>0.05), LD-pruned variants.  |
| supp/sex.mt.chr.variants.txt  | List of 5,574 chrX, chrY and MT variants.  |
| supp/strata.merging.scheme_cohorts.txt | A file showing how the 42 cohorts were merged into 13 strata.  |
