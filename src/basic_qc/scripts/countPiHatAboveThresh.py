#! /usr/bin/python

import sys

genome = sys.argv[1]
thresh = 0.1
X = {}

print "FID", "IID", "n"
n = 0
for line in open(genome, 'r'):
    n += 1
    if n == 1:
        continue
    fields = line.split()
    fid1 = fields[0]
    iid1 = fields[1]
    fid2 = fields[2]
    iid2 = fields[3]

    pihat = float(fields[9])

    if pihat > thresh:
        if (fid1,iid1) in X:
            X[(fid1, iid1)] += 1
        else:
            X[(fid1, iid1)] = 1
            
        if (fid2, iid2) in X:
            X[(fid2, iid2)] += 1
        else:
            X[(fid2, iid2)] = 1


for x in X:
    print x[0], x[1], X[x]
