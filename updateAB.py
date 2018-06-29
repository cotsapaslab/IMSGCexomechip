#! /usr/bin/python

import sys

ab = sys.argv[1]

for line in open(ab, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    x = fields[3]
    x = x.replace("[","")
    x = x.replace("]","")
    y = x.split("/")
    print fields[0], fields[1], fields[2], y[0], y[1]
