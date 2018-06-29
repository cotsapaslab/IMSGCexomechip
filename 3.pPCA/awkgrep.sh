#!/usr/bin/bash

# usage sh awkgrep.sh file1 file2 > out
# file1 has a list of SNPs
# file2 has the data we want to extract. The SNP must be in the first column.

awk '
 # load first file into array indexed by field 1
 NR == FNR {
  file1nr = FNR
  #for (i=1; i<=NF; i++) {
  # file1[$1,i] = $i
  #}
  # store the number of fields for this index
  file1nf[$1] = NF
  next
  }
  {
  file2nr = FNR
  if (!file1nf[$1]) {
   next
  }
  else {
   print $0
  }
 } ' $1 $2
