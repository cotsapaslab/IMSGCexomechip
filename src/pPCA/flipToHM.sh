#!/usr/bin/bash

##### Use this script when there are flipping problems between 2 files that need to be merged. The script compares 2 files with HapMap to define strandness. AT/CG SNPs are tested using 0.4 as MAF threshold.

## Usage:
# sh flipToHM.sh hapmap file out
# hapmap: frequency file from the HapMap. 
# file: frequency file from the dataset.
# out: the name of the output file. 
# NOTE: all files are in 1234 allelic format

## read the options
HM=$1
file=$2
out=$3

## prepare the files to join
# HapMap file
echo "SNP A1 A2 MAF" > $out.__tmp1__
awk 'NR>1 {print $2, $3, $4, $5}' $HM |  sort  >> $out.__tmp1__

# file
echo "SNP A1 A2 MAF" > $out.__tmp2__
awk 'NR>1 {print $2, $3, $4, $5}' $file | sort  >> $out.__tmp2__

## join the 2 temp files, the HM one should be first
perl ~nikolaos/codes/merge_tables.pl $out.__tmp1__ $out.__tmp2__ SNP | grep -w -v "NA"  > $out.__tmp3__

## algorithm to make the call

awk 'NR>1 {
	# check for AT/CG SNPs. These have to add to 5.
	if($2+$3==5){
		# check the MAF for SNPs with MAF<0.4. These wil be flipped if necessary. 
		if($4<0.4 || $7<0.4){
			if($2==$5){
				print $1, "OK";
			} else {
				print $1, "FLIP";
			}
		} else {
			print $1, "R";
		}
	} else {
		if(($2==$5 && $3==$6) || ($2==$6 && $3==$5)) {
		 print $1, "OK"
		} else if (($2==(5-$5) && $3==(5-$6)) || ($2==(5-$6) && $3==(5-$5))) {
		 print $1, "FLIP"
		} else {
		 print $1, "R"
		}
	}
}' $out.__tmp3__ > $out.__tmp4__


## make the output files. 2 files:
# flip file. These SNPs will be flipped
awk '$2=="FLIP" {print $1}' $out.__tmp4__ > $out.flip

# remove file. This one has the AT/CG that have MAF>=0.4. These ones are difficult to flip. Better to remove and impute.
# also this one will have SNPs that have mismatching alleles and the ones that have no frequency information
awk '$2=="R" {print $1}' $out.__tmp4__ > $out.remove
awk '$5=="NA" {print $2}' $file >> $out.remove
awk '$5=="NA" {print $2}' $HM >> $out.remove

# make a fils with SNPs that are OK
awk '$2=="OK" {print $1}' $out.__tmp4__ > $out.ok

# give descriptives
n=`wc -l $out.__tmp3__ | awk '{print $1}'`
f=`wc -l $out.flip | awk '{print $1}'`
ok=`wc -l $out.ok | awk '{print $1}'`
r=`wc -l $out.remove | awk '{print $1}'`
echo "The common SNPs between the files are: "$n
echo "SNPs to flip: "$f
echo "SNPs to remove (AT/CG with MAF>=0.4) and/or SNPs with mismatching alleles: "$r
echo "SNPs that are OK: "$ok


# clean the temp files 
rm $out.__tmp*__
