#!/usr/bin/bash

#usage: sh ~nikolaos/codes/runpca1KG.sh filename n k
#filename is the name of the ped file without the ped extension, e.g. if the file is "input.ped" then filename sohlud be "input"
#n is the number of PCAs smartpca is going to output
#k is the number of individuals to exclude. Set to "0" if you don't want to remove any.
#This script assumes that the everything is merged (dataset and 1KG) in one ped file and the pedind file has "1KG" for all 1KG samples.  
#LAST UPDATE: 2010-06-14 17:28:03

file=$1
n=$2
k=$3

PLINK=/broad/dejagerlab/programs/plink
CONVERTF=/fg/wgas/software/nickp/eig2b/old/bin/convertf
SMARTPCA="/fg/wgas/software/nickp/eig3.0/bin/smartpca"
EIGENSTRAT=/fg/wgas/software/nickp/eig2b/old/bin/eigenstrat

set PATH="/fg/wgas/software/nickp/eig2b/old/bin/:$PATH"

# make the par file for convertf
echo "genotypename:    $file.ped" > $file$n$k.convertf.par
echo "snpname:         $file.map" >> $file$n$k.convertf.par
echo "indivname:       $file.pedind" >> $file$n$k.convertf.par
echo "outputformat:    EIGENSTRAT" >> $file$n$k.convertf.par
echo "genotypeoutname: $file$n$k.geno" >> $file$n$k.convertf.par
echo "snpoutname:      $file$n$k.snp" >> $file$n$k.convertf.par
echo "indivoutname:    $file$n$k.ind" >> $file$n$k.convertf.par
echo "familynames:     NO" >> $file$n$k.convertf.par

#run$file$n$k.convertf
$CONVERTF -p $file$n$k.convertf.par 
			
# make par file for smartpca			
echo "1KG" > 1KG
echo "genotypename:    $file$n$k.geno" > $file$n$k.smartpca.par
echo "snpname:         $file$n$k.snp" >> $file$n$k.smartpca.par
echo "indivname:       $file$n$k.ind" >> $file$n$k.smartpca.par
echo "evecoutname:     $file$n$k.evec" >> $file$n$k.smartpca.par
echo "evaloutname:     $file$n$k.eval" >> $file$n$k.smartpca.par
echo "altnormstyle:    NO" >> $file$n$k.smartpca.par
echo "numoutevec:      $n" >> $file$n$k.smartpca.par
echo "numoutlieriter:  $k"  >> $file$n$k.smartpca.par
echo "poplistname:     1KG" >> $file$n$k.smartpca.par 

$SMARTPCA -p $file$n$k.smartpca.par > $file$n$k.smartpca.log
