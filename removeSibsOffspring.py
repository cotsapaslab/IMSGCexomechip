#! /usr/bin/python

import sys

sibs = sys.argv[1]
offspring = sys.argv[2]

X = [sibs, offspring]
for x in X:
    for line in open(x, 'r'):
        line = line.replace("\n", "")
        fields = line.split()
        print fields[0], fields[1]
        print fields[2], fields[3]
