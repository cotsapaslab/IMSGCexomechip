#! /usr/bin/python

import sys
import os

file = sys.argv[1]
cohorts = sys.argv[2]
pheno = sys.argv[3]

C = {}
if os.path.exists(cohorts) == True:
    for line in open(cohorts, 'r'):
        line = line.replace("\n", "")
        fields = line.split()
        C[(fields[0],fields[1])] = fields[2:]

P = {}
if os.path.exists(pheno) == True:
    for line in open(pheno, 'r'):
        line = line.replace("\n", "")
        fields = line.split()
        P[(fields[0],fields[1])] = fields[2:]

if C == {}:
    a = 0
else:
    a = max([len(r) for r in C.values()])

if P == {}:
    b = 0
else:
    b = max([len(r) for r in P.values()])

for line in open(file, 'r'):
    line = line.replace("\n", "")
    fields = line.split("\t")
    fid = fields[0]
    iid = fields[1]
    if (fid,iid) in C:
        c = C[(fields[0],fields[1])]
    else:
        c = []
        for i in range(a):
            c.append("Unknown")            
    if (fid,iid) in P:
        p = P[(fields[0],fields[1])]
    else:
        p = []
        for i in range(b):
            p.append("Unknown")

    for i in range(len(c)):
        fields.append(c[i])
    for i in range(len(p)):
        fields.append(p[i])        

    print "\t".join(fields)
