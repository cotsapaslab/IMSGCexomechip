#! /usr/bin/python

import sys

fam = sys.argv[1]
cohort = sys.argv[2]

C = {}
for line in open(cohort, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    C[(fields[0],fields[1])] = fields[2]


for line in open(fam, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    if (fields[0],fields[1]) in C:
        c = C[(fields[0],fields[1])]
    else:
        c = "Unknown"
    print fields[0], fields[1], fields[2], fields[3],fields[4], c
