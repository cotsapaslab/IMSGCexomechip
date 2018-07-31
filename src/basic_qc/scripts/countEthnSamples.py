#! /usr/bin/python

import sys

fam = sys.argv[1]
file = sys.argv[2]

F = {}
for line in open(fam, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    F[(fields[0],fields[1])] = ""

n = 0
for line in open(file, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    if (fields[0],fields[1]) in F:
        n += 1

print n
    
