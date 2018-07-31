#! /usr/bin/python
import sys

root = sys.argv[1]
n = int(sys.argv[2])

X = {}
for i in range(1,n + 1):
    k = "x" + str(i)
    X[k] = open(root + "." + str(i) + ".HLA.dosage", 'r')

line = X["x1"].readline().replace("\n", "")
while line != "":
    out = line.split()
    for i in range(2,n + 1):
        q = X["x"+str(i)].readline().replace("\n", "")
        fields = q.split()
        for j in range(3, len(fields)):
            out.append(fields[j])
    print " ".join(out)
    line = X["x1"].readline().replace("\n", "")    
