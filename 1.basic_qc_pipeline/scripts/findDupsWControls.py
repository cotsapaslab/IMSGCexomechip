#! /usr/bin/python

import sys

fam1 = sys.argv[1]
fam2 = sys.argv[2]

X = {}
for line in open(fam2, 'r'):
    fields = line.split()
    X[(fields[0],fields[1])] = ""


for line in open(fam1, 'r'):
    fields = line.split()
    if (fields[0],fields[1]) in X:
        print fields[0],fields[1]
    
