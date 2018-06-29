#! /usr/bin/python

import sys
import random

fam = sys.argv[1]
X = []
dir = sys.argv[2]
j=int(sys.argv[3])

for line in open(fam, 'r'):
    line = line.replace("\n", "")
    fields = line.split()

    fid = fields[0]
    iid = fields[1]

    X.append((fid,iid))

random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)
random.shuffle(X)

n = 0
q = 1
out = open("testfam_"+ str(q) + ".fam",'w')

for i in X:
    if n % j == 0:
        if out != None:
            out.close()
        out = open(dir + "keep_"+ str(q) + ".txt",'w')
        q += 1
    n += 1
    
    out.write(i[0] + " " + i[1] + "\n")
    
