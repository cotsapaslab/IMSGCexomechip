## This script calculates identity-by-missingness matrix for multiple datasets to detect possible correlations in patterns of missing data

# Remove variants on X,Y & MT chrs
ls *.QCed.bed | sed -e 's/.bed$//g' > files.txt
for i in $(cat files.txt) ; do
	./plink --bfile $i \
	--exclude sex.mt.chr.variants.txt \
	--make-bed \
	--out $i'.autosomal'\
        --noweb
done

ls *.QCed.autosomal.bed | sed -e 's/.bed$//g' > files1.txt
for i in $(cat files1.txt) ; do
	./plink --bfile $i \
	--exclude hla.variants.extended.txt \
	--make-bed \
	--out $i'.nohla'\
        --noweb
done

## Calculate the missingness rates
ls *QCed.autosomal.nohla.bed | sed -e 's/.bed$//g' > files2.txt
for i in $(cat files2.txt) ; do
	./plink --allow-no-sex \
	--bfile $i \
	--missing \
	--out $i
    --noweb
done

# Produce the missingness per IND file in tab format 
ls *.QCed.autosomal.nohla.imiss > files3.txt
for i in $(cat files3.txt) ; do
	perl deeplink.pl $i tab
done

## Produce the missingness per SNP file in tab format 
ls *.QCed.autosomal.nohla.lmiss > files4.txt
for i in $(cat files4.txt) ; do
	perl deeplink.pl $i tab
done

## Produce a list of variants with 1 or more missing genotypes and/or count the number of missing sites
ls *.QCed.autosomal.nohla.lmiss.tab > files5.txt
for i in $(cat files5.txt) ; do
#	cat $i | awk '{if ($5 != 0) print $2;}' | wc -l # count only
	cat $i | awk '{if ($5 != 0) print $2;}' > $i'.bad.ibm.variants'
done

# Remove bad variants (common & rare) that are causing IBM/batch effects
ls *.QCed.autosomal.nohla.bed | sed -e 's/.bed$//g' > files6.txt
ls *.QCed.autosomal.nohla.lmiss.tab.bad.ibm.variants > files7.txt
paste files6.txt files7.txt |
while read i j; do
	./plink --bfile $i \
	--exclude $j \
	--make-bed \
	--out $i'.bad.out'\
    --noweb
done

# Cluster the individuals based on missing genotypes (i.e. who's who by proportion of missing genotypes)
ls *.QCed.autosomal.nohla.bad.out.bed | sed -e 's/.bed$//g' > files8.txt
for i in $(cat files8.txt) ; do
 bsub -env arg1=$i -o ibm.prep.out -R "rusage[mem=12]" "./plink --bfile $arg1 --cluster-missing --out $arg1"  #LSF
done	


# Remove the intermediate files
rm files*.txt
rm *.QCed.autosomal.nohla.lmiss
rm *.QCed.autosomal.nohla.imiss
rm  *QCed.autosomal.nohla.bed
rm *.QCed.autosomal.bed
