#! /usr/bin/python

import sys

controls = sys.argv[1]
pheno = sys.argv[2]
fam = sys.argv[3]

nPheno = 0
P = {}
for line in open(pheno, 'r'):
    fields = line.split()
    P[(fields[0],fields[1])] = fields[2:]
    nPheno = len(fields[2:])

C = {}
for line in open(controls, 'r'):
    fields = line.split()
    C[(fields[0],fields[1])] = ""

n = 0
for line in open(fam, 'r'):
    n += 1    
    line = line.replace("\n","")
    fields = line.split()
    out = fields[:2]
    if n == 1:
        out = ["FID","IID"]
        p = P[("FID","IID")]
        for x in p:            
            out.append(x)
        out.append("control")

    elif (fields[0], fields[1]) in P:
        p = P[(fields[0],fields[1])]
        for x in p:
            out.append(x)
        if p[0] == "1":
            out.append("2")
        else:
            out.append("-9")
            
    elif (fields[0],fields[1]) in C:
        for i in range(nPheno):
            out.append("1")
        out.append("1")

    else:
        for i in range(nPheno):
            out.append("-9")
        out.append("-9")

    print " ".join(out)
