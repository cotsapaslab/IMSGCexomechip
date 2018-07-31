# !/bin/sh
# Compare your alleles to the referential allelic map

cohort=$1

awk 'BEGIN{OFS = "\t"};
function complA(allele) {
 if (allele == "A") a2="T";
 else if (allele == "C") a2="G";
 else if (allele == "G") a2="C";
 else if (allele = "T") a2="A";
 else a2="NA";
  return a2}; 

function flipped(A1, A2, A3, A4) {
 if ((A1==A4) && (A2==A3)) return 1;
 else return 0;
}  

function compleflipped(B1, B2, B3, B4) {
 if ((B1==complA(B4)) && (B2==complA(B3))) return 1;
 else return 0;
}  

NR==1 {print $0, "Difference"}; NR>1 {
# check if one of the content has NAs or both of them
  if (($2=="NA" && $3=="NA") || ($4=="NA" && $5=="NA")) {  
  print $0, "NA";
# check whether A1s and A2s are exactly the same
  } else if ($2==$4 && $3==$5) {
  print $0, "same";
# check whether A1s and A2s are complementary, e.g. A/C becomes T/G
  } else if ($2==complA($4) && $3==complA($5)) {
  print $0, "complementary";
# check whether A1s and A2s are flipped, e.g. A/C becomes C/A 
  } else if (flipped($2, $3, $4, $5)==1) { 
  print $0, "flipped";
# check if any of the A1s/A2s have 0s.
  } else if ($2==0 || $3==0 || $4==0 || $5==0) {
  print $0, "zeroed";
# check whether A1s and A2s are complementary and flipped
  } else if (compleflipped($2,$3,$4,$5)==1) {
  print $0, "compleflipped";
# remaining have to be erroneous
  } else {
  print $0, "error";
 }
}' $cohort > $cohort'.alleleCheck'

# Count and print the resulting categories
awk '
/same$/ {
  ++x
}
END {
        print x " alleles are the same"
}' $cohort'.alleleCheck'

awk '
/complementary$/ {
  ++x
}
END {
        print x " alleles are complementary"
}' $cohort'.alleleCheck'

awk '
/flipped$/ {
  ++x
}
END {
        print x " alleles are flipped"
}' $cohort'.alleleCheck'

awk '
/compleflipped$/ {
  ++x
}
END {
        print x " alleles are complementary and flipped"
}' $cohort'.alleleCheck'

awk '
/zeroed$/ {
  ++x
}
END {
        print x " alleles are zeroed"
}' $cohort'.alleleCheck'

awk '
/error$/ {
  ++x
}
END {
        print x " alleles are erroneous"
}' $cohort'.alleleCheck'

