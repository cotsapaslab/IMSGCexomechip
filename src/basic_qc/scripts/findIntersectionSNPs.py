#! /usr/bin/python

import sys

bims = sys.argv[1:]

Q = {}
for b in bims:
    if b.find(".bim") == -1:
        b = b + ".bim"
    for line in open(b, 'r'):
        fields = line.split()
        snp = fields[1]
        if snp in Q:
            Q[snp] += 1
        else:
            Q[snp] = 1

n = len(bims)
for snp in Q:
    if Q[snp] == n:
        print snp
