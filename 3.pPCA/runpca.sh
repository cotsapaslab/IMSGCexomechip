#!/usr/bin/bash

file=$1
n=$2
k=$3
PLINK=/fg/DeJager_Nikos/programs/plink/plink
CONVERTF=/broad//dejagerlab/programs/eigen/convertf
SMARTPCA="/broad//dejagerlab/programs/eigen/smartpca"
EIGENSTRAT=/broad//dejagerlab/programs/eigen//eigenstrat

set PATH="/broad//dejagerlab/programs/eigen/:$PATH"

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
echo "genotypename:    $file$n$k.geno" > $file$n$k.smartpca.par
echo "snpname:         $file$n$k.snp" >> $file$n$k.smartpca.par
echo "indivname:       $file$n$k.ind" >> $file$n$k.smartpca.par
echo "evecoutname:     $file$n$k.evec" >> $file$n$k.smartpca.par
echo "evaloutname:     $file$n$k.eval" >> $file$n$k.smartpca.par
echo "altnormstyle:    NO" >> $file$n$k.smartpca.par
echo "numoutevec:      $n" >> $file$n$k.smartpca.par
echo "numoutlieriter:  $k"  >> $file$n$k.smartpca.par
echo "snpweightoutname: $file$n$k.weights"  >> $file$n$k.smartpca.par

$SMARTPCA -p $file$n$k.smartpca.par > $file$n$k.smartpca.log
