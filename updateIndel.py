#! /usr/bin/python

import sys

bim = sys.argv[1]

for line in open(bim, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    a = fields[4]
    b = fields[5]

    if fields[1].find("IND") == -1:
        continue
    
    if a == "-":
        i = b
    elif b == "-":
        i = a

    print fields[1], "D", "I", "-", i
