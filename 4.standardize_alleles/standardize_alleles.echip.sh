#! /bin/bash 

########################################################################
# This script updates ExomeChip allelic map to match with 1kG alleles
# CopyLeft 2014 Mitja Mitrovic, Cotsapas Lab, Yale School of Medicine
# Dependencies/required files: 	check_alleles.sh						
#								check_alleles.header					
#								allele_map_comp_flipped_to_original.dat	
#								allele_map_comp_to_original.dat			
#								allele_map_flipped_to_original.dat      
#								echip_original_alleles.dat              
#																		
# Usage: standardize_alleles.echip.sh PLINK.bim								
########################################################################

# Define variables
bimfile=$1
cohortname=`basename $bimfile .QCed.ibm.clean.pca.clean.bim`
plinkname=`basename $bimfile .bim`

# Get your and referential alleles side by side (output columns SNP, A1, A2, refA1, refA2,)
awk 'BEGIN {OFS ="\t"};{print $2,$5,$6}' $bimfile | sort > bim.temp
sort echip_original_alleles.dat > echip_original_alleles.dat.temp
awk -F " " 'BEGIN{while(getline<"bim.temp") a[$1]=1 } ; a[$1] ==1 {print $0 } ' echip_original_alleles.dat.temp > $cohortname'.alleles.txt.temp1'
awk -F " " 'BEGIN{while(getline<"echip_original_alleles.dat.temp") a[$1]=1 } ; a[$1] ==1 {print $0 } ' bim.temp > $cohortname'.alleles.txt.temp2'
paste $cohortname'.alleles.txt.temp1' $cohortname'.alleles.txt.temp2' | awk '{print $1, $2, $3, $5, $6}' > $cohortname'.alleles.txt'

# Add the header
echo "$(cat check_alleles.header $cohortname'.alleles.txt')" > $cohortname'.alleles.txt'

## Compare your alleles and with reference map (outcome: same, flipped, complementary, compleflipped, zeroed, erroneous)
sh check_alleles.sh $cohortname'.alleles.txt'

# Get the list of your flipped alleles 
grep "flipped" $cohortname'.alleles.txt.alleleCheck' | awk '{print $1}' > $cohortname'.alleles.txt.alleleCheck.flipped'

# Get the list of your complementary alleles 
grep "complementary" $cohortname'.alleles.txt.alleleCheck' | awk '{print $1}' > $cohortname'.alleles.txt.alleleCheck.complementary'

# Get the list of your complementary+flipped alleles 
grep "compleflipped" $cohortname'.alleles.txt.alleleCheck' | awk '{print $1}' > $cohortname'.alleles.txt.alleleCheck.compleflipped'

# Get the list of your erroneous alleles (should show indels)
grep "error" $cohortname'.alleles.txt.alleleCheck' | awk '{print $1}' > $cohortname'.alleles.txt.alleleCheck.erroneous'

# Create allelic map for flipped alleles (old A1, old A2, ref A1, ref A2)
awk -F " " 'BEGIN{while(getline<"'$cohortname'.alleles.txt.alleleCheck.flipped''") a[$1]=1 } ; a[$1] ==1 {print $0 } ' allele_map_flipped_to_original.dat > $cohortname'.allele_map_flipped_to_original.txt'

# Create allelic map for complementary alleles (old A1, old A2, ref A1, ref A2)
awk -F " " 'BEGIN{while(getline<"'$cohortname'.alleles.txt.alleleCheck.complementary''") a[$1]=1 } ; a[$1] ==1 {print $0 } ' allele_map_comp_to_original.dat > $cohortname'.allele_map_comp_to_original.txt'

# Create allelic map for complementary and flipped alleles (old A1, old A2, ref A1, ref A2)
awk -F " " 'BEGIN{while(getline<"'$cohortname'.alleles.txt.alleleCheck.compleflipped''") a[$1]=1 } ; a[$1] ==1 {print $0 } ' allele_map_comp_flipped_to_original.dat > $cohortname'.allele_map_comp_flipped_to_original.txt'

# Update the flipped alleles
./plink --bfile $plinkname --update-alleles $cohortname'.allele_map_flipped_to_original.txt' --make-bed --out $plinkname'_flipped' --noweb

# Update the complementary alleles
./plink --bfile $plinkname'_flipped' --update-alleles $cohortname'.allele_map_comp_to_original.txt' --make-bed --out $plinkname'_comp' --noweb

# Update the compleflipped alleles
./plink --bfile $plinkname'_comp' --update-alleles $cohortname'.allele_map_comp_flipped_to_original.txt' --make-bed --out $plinkname'_flipped_comp' --noweb

# Remove vars with erroneous alleles 
./plink --bfile $plinkname'_flipped_comp' --exclude $cohortname'.alleles.txt.alleleCheck.erroneous' --make-bed --out $plinkname'.stand_alleles' --noweb


# Clean-up
rm *.temp*
rm $cohortname.*
rm *_flipped.*
rm *_comp.*
rm *_flipped_comp.*
rm *.hh
rm *.nosex