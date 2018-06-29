#! /usr/bin/python

import sys
import os

badsamples = sys.argv[1]
gapSummary = sys.argv[2]

G = {}
if os.path.exists(gapSummary) == True:
    n = 0
    for line in open(gapSummary, 'r'):
        n += 1
        if n < 6:
            continue
        line = line.replace("\n", "")
        fields = line.split("\t")

        fid = fields[0]
        iid = fields[1]
        chip = fields[48].split("-")[1].split("_")[0]
        G[(fid,iid)] = chip

for line in open(badsamples, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    fid = fields[0]
    iid = fields[1]

    if (fid,iid) in G:
        g = G[(fid,iid)]
    else:
        g = "Unknown"

    fields.append(g)
    print "\t".join(fields)

