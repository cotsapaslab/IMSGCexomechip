#! /usr/bin/python

import sys

genomeResult = sys.argv[1]
fam = sys.argv[2]

X = {}
for line in open(fam, 'r'):
    fields = line.split()
    fid = fields[0]
    iid = fields[1]
    mid = fields[2]
    pid = fields[3]
    X[(fid, iid)] = [mid,pid]



for line in open(genomeResult, 'r'):
    if line.find("FID") != -1:
        continue
    fields = line.split()
    fid1 = fields[0]
    iid1 = fields[1]
    fid2 = fields[2]
    iid2 = fields[3]

    x = X[(fid1,iid1)]
    y = X[(fid2,iid2)]
    if fid2 not in x and iid2 not in x and fid1 not in y and iid1 not in y:
        print fid1,iid1,fid2,iid2
