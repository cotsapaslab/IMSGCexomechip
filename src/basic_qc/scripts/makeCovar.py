#! /usr/bin/python

import sys

file = sys.argv[1]

out = ["FID","IID", "C1", "C2","C3","C4","C5","C6","C7","C8","C9","C10"]
print " ".join(out)

n = 0
for line in open(file, 'r'):
    n += 1
    if n == 1:
        continue

    line = line.replace("\n", "")
    fields = line.split()
    id = fields[0].split(":")
    out = [id[0],id[1]]
    for i in range(1,11):
        out.append(fields[i])

    print " ".join(out)
