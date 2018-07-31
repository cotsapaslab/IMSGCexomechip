#! /usr/bin/python

import sys

bim = sys.argv[1]
fhg19 = sys.argv[2]

X = {}
for line in open(fhg19, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    snp = fields[1]
    a = fields[4]
    b = fields[5]

    X[snp] = [a,b]

for line in open(bim, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    snp = fields[1]
    a = fields[4]
    b = fields[5]
    S = [("A","T"),("T","A"),("G","C"),("C","G")]
    if snp in X:
        if (a,b) in S:
            continue
        elif a == "0" or b == "0":
            continue
        elif X[snp] == ["0","0"]:
            continue
        elif b != X[snp][0] and b != X[snp][1] and a != X[snp][0] and a != X[snp][1]:
            print snp, X[snp], a, b
