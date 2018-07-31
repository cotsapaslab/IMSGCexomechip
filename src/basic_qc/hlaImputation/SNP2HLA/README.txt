#################################################################################################################################
#
# SNP2HLA: Imputation of HLA amino acids and classical alleles from SNP genotypes
#
# Author: Sherman Jia (xiaomingjia@gmail.com)
#         + Small modifications by Buhm Han (buhmhan@broadinstitute.org)
#
# USAGE: ./SNP2HLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink {optional: max_memory[mb] window_size}
#
#################################################################################################################################

Thank you for downloading SNP2HLA. To use this package:

1. Download Plink for your platform (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml). Copy the "plink" run file into the current directory (with SNP2HLA.csh).

2. Download Beagle (version 3.0.4) .jar files into the current directory.
- "beagle.jar" from http://faculty.washington.edu/browning/beagle/beagle.html#download
- "beagle2linkage.jar" and "linkage2beagle.jar" from  http://faculty.washington.edu/browning/beagle_utilities/utilities.html
We recommend downloading version 3.0.4 for the compatability issue, even if it is not the newest version. Beagle web page described above includes links for all past-version binaries.

3. Run SNP2HLA with sample data provided (10 samples from British 1958 Birth Cohort, HapMap CEU reference dataset) using the following command:

tcsh: ./SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED plink 2000 1000
bash: ./SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED ./plink 2000 1000

In the above example,
- 1958BC is the SNP genotype plink files (.bed/.bim/.fam),
- HM_CEU_REF is the reference dataset (.bgl.phased/.markers)
- plink is the pointer to the PLINK software
- 2000 is the maximum java heap size (in mb) for imputation using
Beagle (increase as needed)
- 1000 is the marker window sizw that Beagle uses for phasing and imputation


SNP2HLA will also run with default parameters if memory and window size are not provided:
(java mamory = 2Gb, marker window size = 1000)

tcsh: ./SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED plink 
bash: ./SNP2HLA.csh 1958BC HM_CEU_REF 1958BC_IMPUTED ./plink 

----------------------------------------------------------------------------------

Files included in this package:

1. SNP2HLA.csh: Performs imputation (via Beagle) after SNP QC (using PLINK)
2. Merge_tables.pl: Merges files according to indices in a particular column (called by SNP2HLA.csh)
3. ParseDosage.csh: Converts .gprobs (Beagle) file to .dos (PLINK) file
4. HapMap CEU reference dataset (Plink and Beagle formats)
5. Sample SNP dataset of 10 individuals from Britist 1958 Birth Cohort (1958BC.bed/.bim/.fam)

----------------------------------------------------------------------------------

Input files (for SNP2HLA.csh):

1. SNP dataset (.bed/bim/fam PLINK format)
2. Reference dataset (.bgl.phased/.markers Beagle format)
3. Pointer to plink

Output files:

- {OUTPUT}.dosage: PLINK format dosage data (recommended for downstream analysis)
- {OUTPUT}[.bed/.bim/.fam/.ped/.map]: PLINK format best-guess genotype files
- {OUTPUT}.bgl.gprobs: imputation posterior probabilities for SNPs, HLA alleles, and HLA amino acids
- {OUTPUT}.bgl.phased: imputation best-guess genotypes 
- {OUTPUT}.bgl.r2:     imputation predicted r2 with true genotypes

----------------------------------------------------------------------------------

Marker Nomenclature: For binary encodings, P = Present, A = Absent.

1. Classical HLA alleles: HLA_[GENE]_[ALLELE]. 
- HLA_C_0304 = HLA-C:03:04 (four-digit allele)
- HLA_DRB1_07 = HLA-DRB1:07 (two-digit allele)

2. HLA Amino Acids: AA_[GENE]_[AMINO ACID POSITION]_[GENETIC POSITION]_[ALLELE]. 
- AA_A_56_30018678_G = amino acid 56 of HLA-A, genetic position 30018678 (center of codon), allele = G (Gly) of multi-allelic position
- AA_C_291_31345793 = amino acid 291 of HLA-C, genetic position 31345793, bi-allelic (check {OUTPUT}.bim for alleles, P = Present, A = Absent)

3. HLA intragenic SNPS: SNP_[GENE]_[POSITION]_[ALLELE]
- SNP_B_31430319_G = SNP at position 31430319 of HLA-B, allele = G (guanine) of multi-allelic position
- SNP_DRB1_32659974 = SNP at position 32659974 of HLA-DRB1, bi-allelic (check {OUTPUT}.bim for alleles, P = Present, A = Absent)
- SNP_DQB1_32740666_AT = SNP at position 32740666 of HLA-DQB1, alleles = A (adenine) or T (thymine), (check {OUTPUT}.bim for alleles, P = Present, A = Absent)

4. Insertions / deletions: [VARIANT]_[GENE]_[POSITION]_[INSERTION/x=DELETION]
- AA_C_339_31345102_x = deletion at amino acid 339 in HLA-C, genetic position 31345102 (center of codon), (check {OUTPUT}.bim for alleles, P = deletion Present, A = deletion Absent)
- INS_C_295x296_31345779_VLAVLA = insertion between amino acids 295 and 296 of HLA-C, amino acid sequence inserted = VLAVLA, (check {OUTPUT}.bim for alleles, P = insertion Present, A = insertion Absent)
- SNP_DQA1_32717217_x = deletion at genetic position 32717217 of HLA-DQA1, (check {OUTPUT}.bim for alleles, P = deletion Present, A = deletion Absent)

----------------------------------------------------------------------------------

Association testing in PLINK

plink --noweb --dosage OUTPUT.dosage noheader format=1 --fam OUTPUT.fam --logistic --out OUTPUT.assoc

----------------------------------------------------------------------------------
