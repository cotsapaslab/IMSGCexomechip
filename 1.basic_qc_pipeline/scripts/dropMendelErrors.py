#! /usr/bin/python

import sys

fam = sys.argv[1]
fid = sys.argv[2]

F = {}
for line in open(fid, 'r'):
    line = line.replace("\n", "")
    F[line] = ""


for line in open(fam, 'r'):
    fields = line.split()
    if fields[0] in F:
        print fields[0], fields[1]
