#! /usr/bin/python

import sys

x = sys.argv[1]
y = sys.argv[2]
z = sys.argv[3]

FLIP = {}
for line in open(z, 'r'):
    line = line.replace("\n", "")
    FLIP[line] = ""

X = {}
for line in open(x, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    X[fields[0]] = (fields[3],fields[4])

for line in open(y, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    snp = fields[0]
    if snp in X:
        if snp in FLIP and X[snp] == (fields[3],fields[4]):
            print snp
