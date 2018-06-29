#!/usr/bin/bash
# Please acknowledge Nikolaos Patsopoulos when using pPCA scripts

# usage: sh PCA1_1KG_MSchip_exome_hg19.sh dir gwas population

# hard-wired directories and programs
ref="/path/to/1kG/reference/panel/" # Files are expected to be in PLINK (.bed/.bim/.fam) format
plink="/path/to/plink/"


dir=$1
gwas=$2
pop=$3
# the populations are: EUR ASN AFR AMR


# first rename the SNPs to chr:position
awk '{print $2, $1":"$4}' $dir/$gwas.bim > $dir/"$gwas".snp.ref
for i in bed bim fam
do
  cp $dir/$gwas.$i $dir/$gwas.bak.$i
done

$plink --noweb --bfile $dir/$gwas.bak --update-map $dir/$gwas.snp.ref --update-name --recode --out $dir/$gwas
$plink --noweb --file $dir/$gwas --make-bed --out $dir/$gwas 

# make sure that there no morgans listed. This will produce an error later whilw merging with 1KG data
cp $dir/$gwas.bim $dir/$gwas.bim.bak
awk '{print $1, $2, "0", $4, $5, $6}' $dir/$gwas.bim.bak > $dir/$gwas.bim

# make the frq  file
$plink --noweb --bfile $dir/"$gwas"  --allele1234 --freq --out $dir/"$gwas"
sed 's/\(:[0-9].[^ ]*\):.[^ ]*/\1/' $dir/"$gwas".frq > $dir/"$gwas".frq2

# get the SNPs per chromosome from the 1KG data. The issue here is that the SNP names are not the same. We have to use the chromosome and position combination to do the matching and then rename the 1KG based on the chip's IDs. Another issue will be the set of SNPs that are present also in the exome content. 
for i in {1..22}
do
 awk -v chr=$i '$1==chr {print $1":"$4}' $dir/$gwas.bim | awk '!x[$0]++' > $dir/$gwas.snpsfor1KG.$i
 awk '{print $1":"$4, $2}' $ref/chr$i.$pop.bim > $dir/$gwas.1KG.chr$i.$pop.bim.temp
 sh awkgrep.sh $dir/$gwas.snpsfor1KG.$i $dir/$gwas.1KG.chr$i.$pop.bim.temp | awk '{print $2}' > $dir/$gwas.1KG.chr$i.$pop.bim
done
#rm $dir/$gwas.1KG.chr*.$pop.bim.temp

# make one 1KG file locally
for i in {1..22}
do
 $plink --noweb --bfile $ref/chr$i.$pop --extract $dir/$gwas.1KG.chr$i.$pop.bim --make-bed --out $dir/$gwas.1KG.chr$i.$pop
done 

# check for duplicates. If any then remove them
for i in {1..22}
do
 awk 'x[$1, $4]++' $dir/$gwas.1KG.chr$i.$pop.bim > $dir/$gwas.1KG.chr$i.$pop.dups
done 
# also remove the I/D
for i in {1..22}
do
 awk '$5=="I" || $6=="I"' $dir/$gwas.1KG.chr$i.$pop.bim >> $dir/$gwas.1KG.chr$i.$pop.dups
done
for i in {1..22}
do
 if [ -s $dir/$gwas.1KG.chr$i.$pop.dups ]
 then
  $plink --noweb --bfile $dir/$gwas.1KG.chr$i.$pop --exclude $dir/$gwas.1KG.chr$i.$pop.dups --make-bed --out $dir/$gwas.1KG.chr$i.$pop
 fi 
done


# join the chromosomes together
if [ -a $dir/$gwas.$pop.merge ]
then
 rm $dir/$gwas.$pop.merge
fi
# create the merge file
for i in {2..22}
do
 echo "$dir/$gwas.1KG.chr$i.$pop.bed $dir/$gwas.1KG.chr$i.$pop.bim $dir/$gwas.1KG.chr$i.$pop.fam" >> $dir/$gwas.$pop.merge
done

$plink --noweb --bfile $dir/$gwas.1KG.chr1.$pop --merge-list $dir/$gwas.$pop.merge --make-bed --out $dir/$gwas.1KG.$pop


# make the freq file for the 1KG. We need to change the SPN IDs a bit
awk '{print $2, $1":"$4}' $dir/$gwas.1KG.$pop.bim > $dir/$gwas.1KG.$pop.snp_update1
awk '{print $1":"$4, $2}' $dir/$gwas.1KG.$pop.bim > $dir/$gwas.1KG.$pop.snp_update2
$plink --noweb --bfile $dir/$gwas.1KG.$pop --update-map $dir/$gwas.1KG.$pop.snp_update1 --update-name --recode --out $dir/$gwas.1KG.$pop
$plink --noweb --bfile $dir/$gwas.1KG.$pop  --allele1234 --freq --out $dir/$gwas.1KG.$pop
$plink --noweb --file $dir/$gwas.1KG.$pop --make-bed --out $dir/$gwas.1KG.$pop
$plink --noweb --bfile $dir/$gwas.1KG.$pop  --allele1234 --freq --out $dir/$gwas.1KG.$pop


# run the flipping script
sh flipToHM.sh $dir/$gwas.1KG.$pop.frq $dir/"$gwas".frq2 $dir/$gwas.1KG

mv $dir/"$gwas".bim $dir/"$gwas".bim.bak 
sed 's/\(:[0-9].[^ ]*\):.[^\t]*/\1/' $dir/"$gwas".bim.bak > $dir/"$gwas".bim


# flip SNPs
$plink --noweb --bfile $dir/"$gwas" --flip $dir/$gwas.1KG.flip --make-bed --out $dir/$gwas.2
# remove SNPs
$plink --noweb --bfile $dir/$gwas.2 --exclude $dir/$gwas.1KG.remove --alleleACGT --make-bed --out $dir/$gwas.3
#$plink --noweb --bfile $dir/$gwas.3 --update-map $dir/1KG.$pop.snp_update2 --update-name --recode --out $dir/$gwas.3
#$plink --noweb --file $dir/$gwas.3 --make-bed --out $dir/$gwas.3


# keep only the commons SNPs from 1KG.NOTE: for this one we will use the ALL files
for i in {1..22}
do
 awk -v chr=$i '$1==chr {print $1":"$4}' $dir/$gwas.3.bim | awk '!x[$0]++' > $dir/$gwas.3.snpsfor1KG.$i
 awk '{print $1":"$4, $2}' $ref/chr$i.ALL.bim > $dir/$gwas.1KG.chr$i.ALL.bim.temp
 sh awkgrep.sh $dir/$gwas.3.snpsfor1KG.$i $dir/$gwas.1KG.chr$i.ALL.bim.temp | awk '{print $2}' > $dir/$gwas.1KG.chr$i.ALL.bim
 $plink --noweb --bfile $ref/chr$i.ALL --alleleACGT --extract $dir/$gwas.1KG.chr$i.ALL.bim --make-bed --out $dir/$gwas.1KG.$i
done
#rm $dir/$gwas.3.ped $dir/$gwas.3.map



# check for duplicates. If any then remove them
for i in {1..22}
do
 awk 'x[$1, $4]++' $dir/$gwas.1KG.$i.bim > $dir/$gwas.1KG.$i.dups
 # also remove the I/D
 awk '$5=="I" || $6=="I"' $dir/$gwas.1KG.$i.bim >> $dir/$gwas.1KG.$i.dups
 if [ -s $dir/$gwas.1KG.$i.dups ]
 then
  $plink --noweb --bfile $dir/$gwas.1KG.$i --exclude $dir/$gwas.1KG.$i.dups --make-bed --out $dir/$gwas.1KG.$i
 fi
done






# join the chromosomes together
if [ -a $dir/ALL.merge ]
then
 rm $dir/ALL.merge
fi
# create the merge file
for i in {2..22}
 do echo "$dir/$gwas.1KG.$i.bed $dir/$gwas.1KG.$i.bim $dir/$gwas.1KG.$i.fam" >> $dir/ALL.merge
done
$plink --noweb --bfile $dir/$gwas.1KG.1 --merge-list $dir/ALL.merge --make-bed --out $dir/$gwas.1KG







# join
awk '{print $2, $1":"$4}' $dir/$gwas.1KG.bim > $dir/$gwas.1KG.snp_update1
awk '{print $1":"$4, $2}' $dir/$gwas.1KG.bim > $dir/$gwas.1KG.snp_update2
$plink --noweb --bfile $dir/$gwas.1KG --update-map $dir/$gwas.1KG.snp_update1 --update-name --recode --out $dir/$gwas.1KG
$plink --noweb --file $dir/$gwas.1KG --make-bed --out $dir/$gwas.1KG
# one last check that the SNPs are the same
$plink --noweb --bfile $dir/$gwas.1KG  --extract $dir/$gwas.3.bim --make-bed --out $dir/$gwas.1KG
$plink --noweb --bfile $dir/$gwas.3 --extract $dir/$gwas.1KG.bim --make-bed --out $dir/$gwas.3b
$plink --noweb --bfile $dir/$gwas.1KG --bmerge $dir/$gwas.3b.bed $dir/$gwas.3b.bim $dir/$gwas.3b.fam --make-bed --out $dir/$gwas.3.1KG_ALL
#rm $dir/$gwas.1KG.???

# check if the merge was succesful. If not then report it and try to fix the issue. I founf out that extremely few exome SNPs have wrong alleles (maybe the position is wrong?)
if [ -a $dir/$gwas.3.1KG_ALL-merge.missnp ]
then
 # how many
 n_miss=`wc -l $dir/$gwas.3.1KG_ALL-merge.missnp | awk '{print $1}'`
 echo "There are "$n_miss" SNPs mismatches between the input files and the 1KG data. These will be removed to proceed, however you should check the file: $dir/$gwas.3.1KG_ALL-merge.missnp." 
 $plink --noweb --bfile $dir/$gwas.1KG  --exclude $dir/$gwas.3.1KG_ALL-merge.missnp --make-bed --out $dir/$gwas.1KG
 $plink --noweb --bfile $dir/$gwas.3 --exclude $dir/$gwas.3.1KG_ALL-merge.missnp --make-bed --out $dir/$gwas.3b
 $plink --noweb --bfile $dir/$gwas.1KG --bmerge $dir/$gwas.3b.bed $dir/$gwas.3b.bim $dir/$gwas.3b.fam --make-bed --out $dir/$gwas.3.1KG_ALL

fi

# make the ped and pedind files
$plink --noweb --bfile $dir/$gwas.3.1KG_ALL --recode --out $dir/$gwas.3.1KG_ALL
sed -e 's/-9$/1KG/' $dir/$gwas.3.1KG_ALL.fam > $dir/$gwas.3.1KG_ALL.pedind
## run smartpca with 1KG projection
sh runpca1KG.sh $dir/$gwas.3.1KG_ALL 30 0
## make the pops file
echo "FID IID POP TYPE" > $dir/$gwas.3.pops
awk 'NR>1 {print $1, $2, $3, "1KG" }' /broad/dejagerlab/1K/phaseIv3/1KG.pops >> $dir/$gwas.3.pops
awk '{print $1, $2, "DATA", $NF}' $dir/$gwas.3.fam >> $dir/$gwas.3.pops
## plot in R
awk 'NR>1' $dir/$gwas.3.1KG_ALL300.evec > $dir/$gwas.3.1KG_ALL300.evecR
R CMD BATCH -$dir/$gwas.3.1KG_ALL300.evecR -$dir/$gwas.3.pops -$dir/$gwas.3.1KG_ALL300 pcaplot1KGbin.R



## run smartpca to remove individuals
# make the ped and pedind files
$plink --noweb --bfile $dir/$gwas.3b --recode --out $dir/$gwas.3b
cp $dir/$gwas.3b.fam  $dir/$gwas.3b.pedind
## run smartpca with out 1KG
# with no outlier removal
sh runpca.sh $dir/$gwas.3b 30 0
# make the R plots and test which PCAs are statistically significant
awk 'NR>1' $dir/$gwas.3b300.evec > $dir/$gwas.3b300.evecR
R CMD BATCH -$dir/$gwas.3b300.evecR -$dir/$gwas.3b300 pcaplotBin.R &

# run wih 10 iterations to exclude individuals
sh runpca.sh $dir/$gwas.3b 30 10

sleep 20
## descreptives of excluded per iteration
n=0
for i in {1..10}
do 
 echo "Iteration "$i 
 awk -v i=$i '$1=="REMOVED" && $5==i {print $3}'  $dir/$gwas.3b3010.smartpca.log | sort -u  > $dir/$gwas.3b.iteration"$i".removed 
 if [ -s $dir/$gwas.3b.iteration"$i".removed ] 
 then 
  n=$(( $n +1)) 
 fi 
done

## remove the individuals per iteration
awk '{print $2, $1}' $dir/$gwas.3b.fam > $dir/$gwas.3b.fam.tmp
for (( i = 1; i <= $n ; i++ )) 
do
 if [ -s $dir/$gwas.3b.iteration"$i".removed ]
 then
  sh codes/awkgrep.sh $dir/$gwas.3b.iteration"$i".removed $dir/$gwas.3b.fam.tmp  | awk '{print $2, $1}' > $dir/$gwas.3b.iteration"$i".remove
  if [ $i -eq 1 ]  
  then
   $plink --noweb --bfile $dir/$gwas.3b --remove $dir/$gwas.3b.iteration"$i".remove --make-bed --out $dir/$gwas.3b."$i"
  else
   j=$(($i-1))
   $plink --noweb --bfile $dir/$gwas.3b."$j" --remove $dir/$gwas.3b.iteration"$i".remove --make-bed --out $dir/$gwas.3b."$i"
  fi
 fi 
done
rm $dir/$gwas.3b.fam.tmp

## join again with 1KG, run smartpca and plot in R
#$plink --noweb --bfile $dir/$gwas.1KG --update-map $dir/1KG.$pop.snp_update1 --update-name --recode --out $dir/$gwas.1KG
#$plink --noweb --file $dir/$gwas.1KG --make-bed --out $dir/$gwas.1KG
#$plink --noweb --bfile $dir/$gwas.1KG --bmerge $dir/$gwas.3b.bed $dir/$gwas.3b.bim $dir/$gwas.3b.fam --make-bed --out $dir/$gwas.3b.1KG_ALL
#rm $dir/$gwas.1KG.???



for (( i = 1; i <= $n ; i++ )) 
do
  # calculate lambda 
  #$plink --noweb --bfile $dir/$gwas.3."$i" --assoc --adjust --out $dir/$gwas.3.$i
  echo -e "\n-----------------"
  echo "Iteration "$i", "`wc -l $dir/$gwas.3b.iteration"$i".remove |awk '{print $1}'` " removed."
  # remove the individuals from the 1KG-merged files
  if [ $i -eq 1 ]
  then
   $plink --noweb --bfile $dir/$gwas.3.1KG_ALL --remove $dir/$gwas.3b.iteration"$i".remove --make-bed --out $dir/$gwas.3b.1KG_ALL."$i"
   #  make the files
  else
   j=$(($i-1))
   $plink --noweb --bfile $dir/$gwas.3b.1KG_ALL."$j" --remove $dir/$gwas.3b.iteration"$i".remove --make-bed --out $dir/$gwas.3b.1KG_ALL."$i"
  fi 
  $plink --noweb --bfile $dir/$gwas.3b.1KG_ALL."$i" --recode --out $dir/$gwas.3b.1KG_ALL."$i"
  sed -e 's/-9$/1KG/' $dir/$gwas.3b.1KG_ALL."$i".fam > $dir/$gwas.3b.1KG_ALL."$i".pedind
  # run pca I need to bsub this and monitor if it's done
  sh runpca1KG.sh $dir/$gwas.3b.1KG_ALL."$i" 30 0 
done
rm $dir/$gwas.3b.1KG_ALL.*.ped $dir/$gwas.3b.1KG_ALL.*.map


for (( i = 1; i <= $n ; i++ ))
do
 echo -e "\n-----------------"
 echo "Iteration "$i", "`wc -l $dir/$gwas.3b.iteration"$i".removed |awk '{print $1}'` " removed."
done

## make the R plots
for (( i = 1; i <= $n ; i++ )) 
do 
 awk 'NR>1' $dir/$gwas.3b.1KG_ALL."$i"300.evec > $dir/$gwas.3b.1KG_ALL."$i"300.evecR 
 R CMD BATCH -$dir/$gwas.3b.1KG_ALL."$i"300.evecR  -$dir/$gwas.3.pops -$dir/$gwas.3b.1KG_ALL."$i"300 pcaplot1KGbin.R 
done
 
## calculate the PCAs without the 1KG data.
for (( i = 1; i <= $n ; i++ ))
do
 # make the ped and pedind files
 $plink --noweb --bfile $dir/"$gwas".3b.$i --recode --out $dir/"$gwas".3b.$i
 cp $dir/"$gwas".3b.$i.fam  $dir/"$gwas".3b.$i.pedind
 sh runpca.sh $dir/$gwas.3b.$i 30 0
 # make the R plots and test which PCAs are statistically significant
 awk 'NR>1' $dir/$gwas.3b."$i"300.evec > $dir/$gwas.3b."$i"300.evecR
 R CMD BATCH -$dir/$gwas.3b."$i"300.evecR -$dir/$gwas.3b."$i"300 pcaplotBin.R &
done



# restore the bim file
mv $dir/"$gwas".bim.bak $dir/"$gwas".bim



## reminder to evaluate the plots manually
echo "CHECK THE PCA PLOTS TO SEE IF OTHER SAMPLES HAVE TO BE REMOVED"
