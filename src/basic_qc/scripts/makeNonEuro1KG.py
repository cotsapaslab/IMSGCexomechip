#! /usr/bin/python

import sys

cohort = sys.argv[1]

EUR = ["FIN","GBR","CEU","IBS","TSI","IBD","DILI","NIMH_AUT","NIMH_SC"]

for line in open(cohort, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    if fields[2] not in EUR:
        print fields[0], fields[1]

