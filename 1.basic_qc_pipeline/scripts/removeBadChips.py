#! /usr/bin/python

import sys

file = sys.argv[1]
gapSummary = sys.argv[2]

X = {}
for line in open(file, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    chip = fields[2]
    if chip in X:
        X[chip] += 1
    else:
        X[chip] = 1

n = 0
for line in open(gapSummary, 'r'):
    n += 1
    if n < 6:
        continue
    line = line.replace("\n", "")
    fields = line.split("\t")
    chip = fields[48].split("-")[1].split("_")[0]

    if chip in X:
        if X[chip] > 6: # greater than 50% of chip failed
            print fields[0], fields[1]
        
    
