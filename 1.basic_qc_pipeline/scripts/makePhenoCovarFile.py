#! /usr/bin/python

import sys

pheno = sys.argv[1]
covar = sys.argv[2]
fam = sys.argv[3]
p = sys.argv[4]

P = {}
n = 0
for line in open(pheno, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    n += 1
    if n == 1:
        k = fields.index(p)
        continue
    else:
        if fields[k] == "-9":
            fields[k] = "NA"
        elif fields[k] == "1":
            fields[k] = "0"
        elif fields[k] == "2":
            fields[k] = "1"
        P[(fields[0],fields[1])] = fields[k]

C = {}
n = 0
for line in open(covar, 'r'):
    line = line.replace("\n", "")
    n += 1
    fields = line.split()
    if n == 1:
        continue
    else:
        C[(fields[0],fields[1])] = fields[2:]

print "FID IID pheno C1 C2 C3 C4 C5 C6 C7 C8 C9 C10"
for line in open(fam, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    out = [fields[0],fields[1]]
    if (fields[0],fields[1]) in P:
        out.append(P[(fields[0],fields[1])])
    else:
        out.append("NA")

    if (fields[0],fields[1]) in C:
        c = C[(fields[0],fields[1])]
        for i in range(len(c)):
            out.append(c[i])
    else:
        for i in range(10):
            out.append("NA")

    print " ".join(out)
