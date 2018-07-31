#!/bin/bash

## Check for Errors

function die(){
    echo "Error in SNP2HLA.sh" >> $LOG
    exit 1
    }
#trap die ERR

if [ $(($#)) -lt 4 ]; then
    echo "USAGE: ./ImputeHLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink {optional: java_max_memory[mb] marker_window_size}" >> $LOG
    exit 1
fi

SCRIPTPATH="/humgen/atgu1/fs03/wip/jackie/ecPipelineVers4/hlaImputation/SNP2HLA"

MERGE=$SCRIPTPATH/merge_tables.pl
PARSEDOSAGE=$SCRIPTPATH/ParseDosage.csh

# CHECK FOR DEPENDENCIES
if [ ! -e $SCRIPTPATH/beagle.jar ]; then
    echo "Please install Beagle 3 (http://faculty.washington.edu/browning/beagle/beagle.html#download) and copy the run file (beagle.3.0.4/beagle.jar) into $SCRIPTPATH/"
    exit 1
elif [ ! -e $SCRIPTPATH/linkage2beagle.jar ]; then
    echo "Please copy linkage2beagle.jar (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) into $SCRIPTPATH/"
    exit 1
elif [ ! -e $SCRIPTPATH/beagle2linkage.jar ]; then # We use beagle2linkage (Buhm, 8/13/12)
    echo "Please copy beagle2linkage.jar (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/beagle2linkage.jar) into $SCRIPTPATH/"
    exit 1
elif [ ! -e $MERGE ]; then
    echo "Please copy merge_tables.pl (included with this package) into $SCRIPTPATH/"
    exit 1
elif [ ! -e $PARSEDOSAGE ]; then
    echo "Please copy ParseDosage.csh (included with this package) into $SCRIPTPATH/"
    exit 1
fi

# INPUTS
INPUT=$1
REFERENCE=$2
OUTPUT=$3
KEEP=$4

MEM=8000 # Default java memory 2000 Mb (2Gb)
WIN=1000 # Default Beagle phasing/imputation window = 1000 markers

JAVATMP=$OUTPUT.javatmpdir
mkdir -p $JAVATMP
alias plink='plink --noweb --silent --allow-no-sex'
alias beagle='java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/beagle.jar'
alias linkage2beagle='java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/linkage2beagle.jar'
alias beagle2linkage='java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/beagle2linkage.jar'

# Functions to run
EXTRACT_MHC=1
FLIP=1
CONVERT_IN=1
IMPUTE=1
CONVERT_OUT=1
CLEANUP=0

# SET PARAMETERS
TOLERATED_DIFF=.15
i=1

MHC=$OUTPUT.MHC

if [ -n $KEEP ]; then
    plink --bfile $INPUT --keep $KEEP --make-bed --out $OUTPUT.tmpKeep
    INPUT=$OUTPUT.tmpKeep
fi

if [ $EXTRACT_MHC -eq 1 ]; then
    echo "[$i] Extracting SNPs from the MHC." >> $LOG; let i++
    plink --bfile $INPUT --chr 6 --from-mb 29 --to-mb 34 --maf 0.025 --make-bed --out $OUTPUT.MHC
fi

if [ $FLIP -eq 1 ]; then
    echo "[$i] Performing SNP quality control." >> $LOG; let i++

    # Identifying non-A/T non-C/G SNPs to flip
    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.bim >> $OUTPUT.tmp1
    echo "SNP 	POSR	A1R	A2R" > $OUTPUT.tmp2
    cut -f2,4- $REFERENCE.bim >> $OUTPUT.tmp2
    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP |  grep -v -w NA > $OUTPUT.SNPS.alleles

    awk '{if ($3 != $6 && $3 != $7){print $1}}' $OUTPUT.SNPS.alleles > $OUTPUT.SNPS.toflip1
    plink --bfile $MHC --flip $OUTPUT.SNPS.toflip1 --make-bed --out $MHC.FLP

    # Calculating allele frequencies
    plink --bfile $MHC.FLP --freq --out $MHC.FLP.FRQ
    sed 's/A1/A1I/g' $MHC.FLP.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp

    mv $OUTPUT.tmp $MHC.FLP.FRQ
    $MERGE $REFERENCE.FRQ.frq $MHC.FLP.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.frq
    sed 's/ /\t/g' $OUTPUT.SNPS.frq | awk '{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}' > $OUTPUT.SNPS.frq.parsed
    
    # Finding A/T and C/G SNPs
    awk '{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}' $OUTPUT.SNPS.frq.parsed > $OUTPUT.SNPS.ATCG.frq

    # Identifying A/T and C/G SNPs to flip or remove
    awk '{if ($10 < $9 && $10 < .15){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toflip2
    awk '{if ($4 > 0.4){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toremove

    # Identifying non A/T and non C/G SNPs to remove
    awk '{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > '$TOLERATED_DIFF'){print $1}}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove

    # Making QCd SNP file
    plink --bfile $MHC.FLP --geno 0.2 --exclude $OUTPUT.SNPS.toremove --flip $OUTPUT.SNPS.toflip2 --make-bed --out $MHC.QC
    plink --bfile $MHC.QC --freq --out $MHC.QC.FRQ
    sed 's/A1/A1I/g' $MHC.QC.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp
    mv $OUTPUT.tmp $MHC.QC.FRQ.frq
    $MERGE $REFERENCE.FRQ.frq $MHC.QC.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.QC.frq

    cut -f2 $OUTPUT.SNPS.QC.frq | awk '{if (NR > 1){print $1}}' > $OUTPUT.SNPS.toinclude

    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.QC.bim >> $OUTPUT.tmp1

    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP | awk '{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}' > $MHC.QC.bim

    # Recoding QC'd file as ped
    plink --bfile $MHC.QC --extract $OUTPUT.SNPS.toinclude --recode --out $MHC.QC

        # Making SNP files
    awk '{print "M " $2}' $MHC.QC.map > $MHC.QC.dat
    cut -f2 $MHC.QC.map > $MHC.snps
    cut -d ' ' -f1-5,7- $MHC.QC.ped > $MHC.QC.nopheno.ped

    # Remove temporary files
    rm $OUTPUT.tmp1 $OUTPUT.tmp2
    rm $MHC.FLP*
    rm $MHC.QC.ped $MHC.QC.map
    rm $OUTPUT.SNPS.*
fi

if [ $CONVERT_IN -eq 1 ]; then
    echo "[$i] Converting data to beagle format." >> $LOG; let i++
    java -jar ${SCRIPTPATH}/linkage2beagle.jar $MHC.QC.dat $MHC.QC.nopheno.ped  > $MHC.QC.bgl
    rm $MHC.QC.nopheno.ped
fi

if [ $IMPUTE -eq 1 ]; then
    echo "[$i] Performing HLA imputation (see $OUTPUT.bgl.log for progress)." >> $LOG; let i++
    java -jar ${SCRIPTPATH}/beagle.jar markers=$REFERENCE.markers unphased=$MHC.QC.bgl phased=$REFERENCE.bgl.phased gprobs=true niterations=10 nsamples=4 missing=0 verbose=true maxwindow=$WIN out=$OUTPUT.IMPUTED log=$OUTPUT.phasing > $OUTPUT.bgl.log
fi

if [ $CONVERT_OUT -eq 1 ]; then
    echo "[$i] Converting posterior probabilities to PLINK dosage format." >> $LOG; let i++
    x=`echo $OUTPUT | awk 'BEGIN {FS="/"};{print $NF}'`
    # Converting .gprobs to .dosage format
    mv $OUTPUT.IMPUTED.$x.MHC.QC.bgl.phased $OUTPUT.bgl.phased
    mv $OUTPUT.IMPUTED.$x.MHC.QC.bgl.gprobs $OUTPUT.bgl.gprobs
    mv $OUTPUT.IMPUTED.$x.MHC.QC.bgl.r2 $OUTPUT.bgl.r2

    $PARSEDOSAGE $OUTPUT.bgl.gprobs > $OUTPUT.dosage

    echo "[$i] Converting imputation genotypes to PLINK .ped format." >> $LOG; let i++
    cat $OUTPUT.bgl.phased | beagle2linkage $OUTPUT.tmp # Buhm
    cut -d ' ' -f6- $OUTPUT.tmp.ped > $OUTPUT.tmp       # Buhm
    paste -d ' ' $MHC.fam $OUTPUT.tmp | tr -d "\015" > $OUTPUT.ped
    cut -f1-4 $REFERENCE.bim > $OUTPUT.map
    cp $MHC.fam $OUTPUT.fam

    # Create PLINK bed file
    plink --ped $OUTPUT.ped --map $OUTPUT.map --make-bed --out $OUTPUT
fi

if [ $CLEANUP -eq 1 ]; then
    rm $OUTPUT.*.bgl.phased
    rm $OUTPUT.*.bgl.r2
    rm $OUTPUT.*.bgl.gprobs
    rm $OUTPUT.*.dosage
    rm $OUTPUT.MHC.*
    rm $OUTPUT.tmp*
    rm $OUTPUT.IMPUTED.*.bgl.phased.phased
    rm $OUTPUT.phasing.log
    rm $OUTPUT.ped
    rm $OUTPUT.map
    rm $OUTPUT.log
    rm -r $JAVATMP
    rm -f plink.log
    echo "DONE!" >> $LOG
    echo "" >> $LOG
fi