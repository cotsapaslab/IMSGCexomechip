#! /usr/bin/python

import sys

ped = sys.argv[1]


out = ["FID", "IID", "nHET", "percHET", "nMISS", "percMISS"]
print "\t".join(out)

for line in open(ped, 'r'):
    fields = line.split()
    fid = fields[0]
    iid = fields[1]
    nHet = 0
    nTotal = 0
    nMiss = 0
    for i in range(6, len(fields), 2):
        if fields[i] == "0" and fields[i + 1] == "0":
            nMiss += 1
        elif (fields[i] != fields[i+1]):
            nHet += 1
        nTotal += 1

    out = [fid, iid, nHet, (float(nHet)/ float(nTotal))*100, nMiss, (float(nMiss) / float(nTotal))*100]
    out = [str(o) for o in out]
    print "\t".join(out)

