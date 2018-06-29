#! /usr/bin/python

import sys

assoc = sys.argv[1]
fisher = sys.argv[2]
logistic = sys.argv[3]
hla = sys.argv[4]
fisherPerm1000 = sys.argv[5]
fisherPermMillion = sys.argv[6]
hlaAssoc = sys.argv[7]
root = sys.argv[8]
evs = "/humgen/atgu1/fs03/wip/elim/evs/evs_counts.txt"

EVS = {}
for line in open(evs, 'r'):
    line = line.replace("\n", "")
    fields = line.split("\t")
    chr = fields[0]
    pos = fields[1]
    gene = fields[2]
    function = fields[3]
    pph = fields[6]
    maf = fields[9]

    EVS[(chr,pos)] = [gene,function,pph,maf]

A = {}
for line in open(assoc, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    snp = fields[1]
    fA = fields[4]
    fU = fields[5]
    A[snp] = [fA,fU]

for line in open(hlaAssoc, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    snp = fields[1]
    fA = fields[4]
    fU = fields[5]
    if snp not in A:
        A[snp] = [fA,fU]

P = {}
for line in open(fisherPerm1000, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    snp = fields[1]
    p = fields[2]
    P[snp] = [p, "1000"]

for line in open(fisherPermMillion, 'r'):
    line = line.replace("\n", "")
    fields = line.split()
    snp = fields[1]
    p = fields[2]
    P[snp] = [p, "1000000"]
    
out = open(root + "_fisher.csv", 'w')
n = 0
for line in open(fisher, 'r'):
    n += 1
    if n == 1:
        out.write(" ".join(["SNP","CHR","BP","F_A","F_U","P","OR","P_PERM", "N_PERM","GENE","FUNCTION","PPH","MAF_EVS"]) + "\n")
        continue
    line = line.replace("\n", "")
    fields = line.split()
    chr = fields[0]
    pos = fields[2]
    snp = fields[1]
    OR = fields[8]
    p = fields[7]
    fA = fields[4]
    fU = fields[5]
    pPerm = P[snp][0]
    nPerm = P[snp][1]
    if fA == "0" and fU == "0":
        continue

    if (chr,pos) in EVS:
        e = EVS[(chr,pos)]
        gene = e[0]
        function = e[1].split(",")[0]
        pph = e[2].split(",")[0]
        maf = e[3]

        if function  == "missense":
            function = "Missense"
        elif function == "coding-synonymous":
            function = "Synonymous"
        elif function == "missense-near-splice":
            function = "Splice"
        elif function == "stop-gained":
            function = "Nonsense"

        if pph == "benign":
            pph = "Benign"
        elif pph == "probably-damaging":
            pph = "Probably-Damaging"
        elif pph == "possibly-damaging":
            pph = "Possibly-Damaging"
    else:
        gene = ""
        function = ""
        pph = ""
        maf = ""

    out.write(" ".join([snp, chr, pos, fA, fU, p, OR, pPerm, nPerm, gene, function, pph, maf]) + "\n")

out.close()

out = open(root + "_logistic.csv", 'w')
n = 0
for line in open(logistic, 'r'):
    n += 1
    if n == 1:
        out.write(" ".join(["SNP","CHR","BP","F_A","F_U","P","OR","GENE","FUNCTION","PPH","MAF_EVS"]) + "\n")
        continue
    line = line.replace("\n", "")
    fields = line.split()
    chr = fields[0]
    pos = fields[2]
    snp = fields[1]
        
    test = fields[4]
    OR = fields[6]
    p = fields[8]

    if p == "NA":
        continue
        
    if test == "ADD":
        fA = A[snp][0]
        fU = A[snp][1]
            
        if (chr,pos) in EVS:
            e = EVS[(chr,pos)]
            gene = e[0]
            function = e[1].split(",")[0]
            pph = e[2].split(",")[0]
            maf = e[3]

            if function  == "missense":
                function = "Missense"
            elif function == "coding-synonymous":
                function = "Synonymous"
            elif function == "missense-near-splice":
                function = "Splice"
            elif function == "stop-gained":
                function = "Nonsense"

            if pph == "benign":
                pph = "Benign"
            elif pph == "probably-damaging":
                pph = "Probably-Damging"
            elif pph == "possibly-damaging":
                pph = "Possibly-Damaging"

        else:
            gene = ""
            function = ""
            pph = ""
            maf = ""
            
        out.write(" ".join([snp, chr, pos, fA, fU, p, OR, gene, function, pph, maf]) + "\n")
out.close()

out = open(root + "_hla.csv", 'w')
n = 0
for line in open(hla, 'r'):
    n += 1
    if n == 1:
        out.write(" ".join(["SNP","CHR","BP","F_A","F_U","P","OR","GENE","FUNCTION","PPH","MAF_EVS"]) + "\n")
        continue
    line = line.replace("\n", "")
    fields = line.split()
    chr = fields[0]
    pos = fields[2]
    snp = fields[1]
        
    test = fields[4]
    OR = fields[6]
    p = fields[8]

    if p == "NA":
        continue
        
    if test == "ADD":
        fA = A[snp][0]
        fU = A[snp][1]
            
        if (chr,pos) in EVS:
            e = EVS[(chr,pos)]
            gene = e[0]
            function = e[1].split(",")[0]
            pph = e[2].split(",")[0]
            maf = e[3]

            if function  == "missense":
                function = "Missense"
            elif function == "coding-synonymous":
                function = "Synonymous"
            elif function == "missense-near-splice":
                function = "Splice"
            elif function == "stop-gained":
                function = "Nonsense"

            if pph == "benign":
                pph = "Benign"
            elif pph == "probably-damaging":
                pph = "Probably-Damging"
            elif pph == "possibly-damaging":
                pph = "Possibly-Damaging"

        else:
            gene = ""
            function = ""
            pph = ""
            maf = ""
            
        out.write(" ".join([snp, chr, pos, fA, fU, p, OR, gene, function, pph, maf]) + "\n")

    
