#! /bin/bash

## Test for different genotype call rates between cases and controls

ls *.bed | sed -e 's/.bed$//g' > files.txt

for i in $(cat files.txt) ; do
	plink --bfile $i \
	--test-missing \
	--out $i
done

## Produce a list of with significantly different vars (p<0.0001)
ls *.missing > files.txt
for i in $(cat files.txt) ; do
	perl run-diffmiss-qc.pl $i
done

## Count "bad" vars

ls *.fail-diffmiss-qc.txt > files.txt
for i in $(cat files.txt) ; do
	echo $i | awk '{split($0,a,"."); print a[1]}'
	cat $i | wc -l
done

rm files.txt