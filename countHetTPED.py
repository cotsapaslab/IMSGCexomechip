#! /usr/bin/python

import sys

tped = sys.argv[1]


for line in open(tped, 'r'):
    n = 0
    line = line.replace("\n", "")
    fields = line.split()
    for i in range(4,len(fields),2):
        if fields[i] != fields[i+1] and fields[i] != "0" and fields[i+1] != "0":
            n += 1


    print fields[1], n
