#! /usr/bin/python

import sys

fam = sys.argv[1]

for line in open(fam, 'r'):
    fields = line.split()
    fid = fields[0]
    iid = fields[1]

    if fid[:2] == "NA":
        if fid.find("_") != -1:
            print fid, iid
            continue

    if iid[:2] == "NA":
        if iid.find("_") != -1:
            print fid, iid
            continue    
