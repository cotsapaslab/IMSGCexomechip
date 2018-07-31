#! /usr/bin/python

import sys
import os

## Need to put in if sex or mendel errors were removed
def extractFailMiss(log):
    n = "NA"
    if os.path.exists(log) != False:        
        for line in open(log, 'r'):
            if line.find("SNPs failed missingness test") != -1:
                fields = line.split()
                n = fields[0]
    return n

def extractFailHWE(log):
    n = "NA"
    if os.path.exists(log) != False:            
        for line in open(log, 'r'):
            if line.find("markers to be excluded based on HWE test") != -1:
                fields = line.split()
                n = fields[0]
    return n

def countMaleFemale(fam):
    x = {"-9":0,"2":0,"1":0,"0":0}
    for line in open(fam, 'r'):
        fields = line.split()
        x[fields[4]] += 1

    return (str(x["1"]),str(x["2"]),str(x["-9"] + x["0"]))
        
def countChr(bim):
    apa = 0
    x = 0
    y = 0
    mt = 0
    for line in open(bim, 'r'):
        fields = line.split()
        chr = int(fields[0])
        if chr < 23 or chr == 25:
            apa += 1
        elif chr == 23:
            x += 1
        elif chr == 24:
            y += 1
        elif chr == 26:
            mt += 1

    return (str(apa), str(x),str(y),str(mt))
        
        
def makeTable(l,outFile):
    for q in l:
        outFile.write(",".join(q) + "\n")
    return

def countndefined(f):
    n = 0
    X = {}
    for line in open(f, 'r'):
        fields = line.split()
        if fields[2] != "0" or fields[3] != "0":
            X[fields[0]] = ""

    return str(len(X))

def wc(f):
    n = 0
    if os.path.exists(f) == False:
        return str(0)
    
    for line in open(f, 'r'):
        n += 1
    return str(n)

def uniq(l):
    x = {}
    for q in l:
        if os.path.exists(q) == False:
            continue
        for line in open(q):
            line = line.replace("\n", "")
            fields = line.split()
            if fields[0] not in x:
                x[fields[0]] = ""
    return str(len(x))

def countOverlap(f,g):
    X = {}
    if os.path.exists(f) == False:
        return str(0)
    
    for line in open(f, 'r'):
        fields = line.split()
        X[fields[0]] = ""
    Y = {}
    for line in open(g, 'r'):
        fields = line.split()
        if fields[0] in X:
            Y[fields[0]] = ""
    return str(len(Y))

def countnRemoved(q, r):
    X = {}
    for i in q:
        if os.path.exists(i) == False:
            continue
        for line in open(i, 'r'):
            line = line.replace("\n", "")
            fields = line.split()
            X[fields[0]] = ""
            X[fields[2]] = ""

    Y = {}
    for line in open(r, 'r'):
        fields = line.split()
        if fields[0] in X:
            Y[fields[0]] = ""

    return str(len(Y))
        
def writeSuppMethods(inputVar, batches):
    OUTPUT_DIR = inputVar["Output"]

    outFile = open(OUTPUT_DIR + "suppMethods.txt",'w')
    outFile.write( "Supplementary Methods for Exome Chip QC Pipeline\n")
    outFile.write( "Please acknowledge Jackie Goldstein when using the pipeline.\n")
    outFile.write("\n")
    outFile.write( "---------------------------------------------------\n")
    outFile.write("\n")

    o = ["Each batch of data with GenCall and zCall file inputs ("]
    q = []
    for b in batches:
        if batches[b]["gencall"] != "":
            q.append(b)
    o.append(",".join(q))
    o.append(") was QC'd separately as described below:")
    outFile.write( "".join(o) + "\n")
    outFile.write("\n")
    o = ["First, samples were QC'd for the percentage of no calls for common SNPs (> 5% MAF) and the percentage of het calls for rare SNPs (< 5% MAF). "]
    q = []
    r = []
    for b in batches:
        if batches[b]["use1KGdef"] == "1":
            q.append(b)
        else:
            r.append(b)
    if len(q) != 0:
        o.append("For batches ")
        o.append(", ".join(q))
        o.append(", MAFs were calculated from a QC'd 1000 genomes dataset.")
    if len(r) != 0:
        o.append("For batches ")
        o.append(", ".join(r))
        o.append(", MAFs were calculated from SNPs that had a call rate greater than 95%.")
        

    o.append(" A table summarizing the number of samples eliminated from each batch and the thresholds used is included below (Table 1).")
    outFile.write( "".join(o) + "\n")
    outFile.write("\n")
    
    outFile.write( "Second, PLINK was used to generate a genome file from which the degree of relatedness for each sample pair was estimated. A table summarizing the number of related samples and the criteria used is in Table 2.\n")
    outFile.write("\n")

    outFile.write( "Third, Sex Checks and Mendel Checks were performed using PLINK. See Table 3 for details.\n")
    outFile.write("\n")

    outFile.write( "Fourth, EIGENSTRAT was used to automatically detect population ancestry (European, African, Asian, Other) for each sample by using 1000 Genomes samples plus the shared controls as population loadings. A table is provided below summarizing the number of samples for each batch (Table 4).\n")
    outFile.write("\n")

    outFile.write( "Finally, SNP QC was done as follows: The autosomal and pseudoautosomal SNPs were QC'd by first removing SNPs with low call rate based on the original Illumina GenCall calls. Thresholds used for each batch are shown in Table 5. Of the SNPs that had a passing GenCall call rate, the zCall call rate per SNP had to be above a certain threshold (Table 5)  and have a HWE p-value less than 1e-5 in the previously defined European subset. The X chromosome SNPs were filtered in the same way as autosomal SNPs, except only females were used for filtering. Y chromosome SNPs were kept if a maximum of 2 males were called as heteroyzgotes. The mitochondrial SNPs were not filtered at all. A summary of the number of SNPs that passed each filtering step by batch and the thresholds used is in Table 5.\n")
    outFile.write("\n")

    q = []
    for b in batches:
        if batches[b]["qc"] != "" and batches[b]["skip"] == "0":
            q.append(b)
            
    o = ["Batches with data that had been previously QC'd ("]
    o.append(", ".join(q))
    o.append(") was left untouched except EIGENSTRAT was used to identify the population ancestry of each sample (Table 3).")

    outFile.write( "".join(o) + "\n")
    outFile.write("\n")

    if inputVar["qc_only"] != "1":
        outFile.write( "The QC'd data and a shared controls set (details in Table 6) consisting of nonrelated samples were then merged using PLINK. This was followed by PCA on the European samples using Eigenstrat with the shared controls as population loadings (Figure 1) and a final check for duplicate samples among batches using PLINK (pi_hat > 0.9) (Table 7). SNPs were only included in the final merged dataset if they passed QC in all datasets including the shared controls set and had a HWE P > 1e-5 across all European samples. The shared controls set had " + wc(inputVar["Script"] + "controls/controls.bim") + " SNPs that had been QC'd using an earlier version of the pipeline.\n")

    outFile.write("\n")
    if inputVar["qc_only"] != "1" and inputVar["skip_hla"] != "1":
        outFile.write( "HLA Imputation was done on the final merged dataset using the SNP2HLA tool (paper will be out soon in PLOS One).\n")
    
    outFile.write("\n")
    outFile.write( "-------------------------------------------------\n")
    outFile.write("\n")
    
    outFile.write( "Table 1: Sample QC Summary\n")
    t = [["Batch", "n","% Het Call Range","Minimum Call Rate","Used 1KG SNP MAFs","# Samples Excluded for Het Calls", "# Samples Excluded for Call Rate","n Total Removed", "n After Sample QC"]]
    for b in batches:
        p = batches[b]
        batchdir = inputVar["Output"] + b + "/"
        if p["gencall"] != "" and p["skip"] == "0":
            q = [b, wc(p["gencall"] + ".fam"), p["minHet"], p["maxHet"], p["maxNoCall"], p["use1KGdef"], wc(batchdir + "bad.rare.txt"), wc(batchdir + "bad.common.txt"), uniq([batchdir + "bad.rare.txt",batchdir + "bad.common.txt"]), str(int(wc(p["gencall"] + ".fam")) - int(uniq([batchdir + "bad.rare.txt",batchdir + "bad.common.txt"])))]            
            t.append(q)
    makeTable(t,outFile)
    outFile.write("\n")
    outFile.write("\n")
    
    outFile.write( "Table 2: Relatedness Check Summary\n")
    t = [["Batch", "nDuplicate Pairs", "n Sibling Pairs", "n Parent Offspring Pairs", "nDups Removed", "nSibs Removed", "nParentOffspring Removed", "Total Initial", "Total Removed", "Total Remaining"]]
    for b in batches:
        p = batches[b]
        batchdir = inputVar["Output"] + b + "/"
        if p["gencall"] != "" and p["skip"] == "0":
            q = [b, wc(batchdir + "dups.txt"), wc(batchdir + "sibs.txt"), wc(batchdir + "offspring.txt"),countnRemoved([batchdir+"dups.txt"],batchdir + "bad.samples.txt"), countnRemoved([batchdir+"sibs.txt"],batchdir + "bad.samples.txt"), countnRemoved([batchdir+"offspring.txt"],batchdir + "bad.samples.txt"), countnRemoved([batchdir+"dups.txt",batchdir+"sibs.txt",batchdir+"offspring.txt"],batchdir+"bad.samples.txt"), str(int(wc(p["gencall"] + ".fam")) - int(uniq([batchdir + "bad.rare.txt",batchdir + "bad.common.txt"])) - int(uniq([batchdir + "sibs.txt",batchdir + "dups.txt",batchdir + "offspring.txt"])))]
            t.append(q)
    
    makeTable(t,outFile)
    outFile.write("\n")
    outFile.write("\n")

    outFile.write( "Table 3: Sex Check and Mendel Error Summary\n")
    t = [["Batch","nSamples", "nSexCheck Errors", "nSamples Removed","nFamilies Defined", "maxMendelErrors Allowed", "nFailing Families", "nSamples Removed", "nSamples Remaining"]]
    for b in batches:
        p = batches[b]        
        if p["skip"] == "0":
            if p["gencall"] != "":
                root = p["gencall"]
                q = [b,str(int(wc(p["gencall"] + ".fam")) - int(uniq([batchdir + "bad.rare.txt",batchdir + "bad.common.txt"])) - int(uniq([batchdir + "sibs.txt",batchdir + "dups.txt",batchdir + "offspring.txt"]))),wc(batchdir+"sexCheckErrors.txt"), countnRemoved([batchdir+"sexCheckErrors.txt"],batchdir + "bad.samples.txt"),countndefined(root + ".fam"),p["maxMendelError"],wc(batchdir + "familyIDs.todrop.txt"), countnRemoved([batchdir + "mendelErrors.txt"], batchdir + "bad.samples.txt"),str(int(wc(p["gencall"] + ".fam")) - int(uniq([batchdir + "bad.samples.txt"])))]
                t.append(q)
    makeTable(t,outFile)
    outFile.write("\n")
    outFile.write("\n")


    outFile.write( "Table 4: PCA Summary\n")
    t = [["Batch", "nEuro", "nAfr", "nAsian", "nOther", "nTotal"]]
    for b in batches:
        p = batches[b]
        batchdir = inputVar["Output"] + b + "/"
        if p["skip"] == "0":
            if p["gencall"] != "":
                root = p["gencall"]
                q = [b, countOverlap(root + ".fam", batchdir + "european.txt"), countOverlap(root + ".fam", batchdir + "african.txt"), countOverlap(root + ".fam", batchdir + "asian.txt"), str(int(wc(root + ".fam")) - int(uniq([batchdir + "bad.samples.txt"])) - int(countOverlap(root + ".fam", batchdir + "european.txt")) - int(countOverlap(root + ".fam", batchdir + "african.txt")) - int(countOverlap(root + ".fam", batchdir + "asian.txt"))),str(int(wc(root + ".fam")) - int(uniq([batchdir + "bad.samples.txt"])))]
            
            elif p["qc"] != "":
                root = p["qc"]            
                q = [b, countOverlap(root + ".fam", batchdir + "european.txt"), countOverlap(root + ".fam", batchdir + "african.txt"), countOverlap(root + ".fam", batchdir + "asian.txt"), str(int(wc(root + ".fam")) - int(countOverlap(root + ".fam", batchdir + "european.txt")) - int(countOverlap(root + ".fam", batchdir + "african.txt")) - int(countOverlap(root + ".fam", batchdir + "asian.txt"))),wc(root + ".fam")]
            t.append(q)
    
    makeTable(t,outFile)
    outFile.write("\n")
    outFile.write("\n")


    outFile.write( "Table 5: SNP QC Summary\n")
    t = [["Batch", "nP+PA initial", "nX initial", "nY initial", "nMT initial", "GenCall Missing Rate", "nSNPs Removed due to GenCall Missing Rate (A+PA)", "zCall Missing Rate", "nSNPs Removed due to zCall Missing Rate (A+PA)", "nSNPs removed due to HWE (A+PA)", "nSNPs Removed due to GenCall Missing Rate (X)", "nSNPs Removed due to zCall Missing Rate (X)", "nSNPs Removed due to HWE (X)", "nSNPs Removed due to n Het Calls (Y)", "nA+PA after QC", "nX after QC", "nY after QC", "nMT", "nTotal"]]
    for b in batches:
        p = batches[b]
        batchdir = inputVar["Output"] + b + "/"
        if p["skip"] == "0":
            if p["gencall"] != "":
                root = p["gencall"]
                apa, x, y, mt = countChr(root + ".bim")
                apaqc, xqc, yqc,mtqc = countChr(batchdir + b + ".bim")

                q = [b, apa, x, y, mt, p["gencall_missrate"], str(int(apa) - int(wc(batchdir + "keepSNPS.gcAuto.txt"))),p["zcall_missrate"], str(int(apa) - int(wc(batchdir + "keepSNPS.gcAuto.txt")) - int(extractFailMiss(batchdir+"zCall.autosomal.eur.log"))), extractFailHWE(batchdir+"zCall.autosomal.eur.log"),str(int(x) - int(wc(batchdir + "keepSNPS.gcX.txt"))), str(int(x) - int(wc(batchdir + "keepSNPS.gcX.txt")) - int(extractFailMiss(batchdir+"zcall.tmp3.log"))), extractFailHWE(batchdir+"zcall.tmp3.log"), str(int(y) - int(wc(batchdir + "keepY.txt"))), apaqc, xqc, yqc, mtqc, wc(batchdir+b+ ".bim")]
                t.append(q)
    
    makeTable(t,outFile)
    outFile.write("\n")
    outFile.write("\n")

    outFile.write( "Table 6: Shared Controls Summary\n")
    t = [["Source", "nEur", "nAfr", "nAsian", "nOther", "nTotal", "Chip Version"]]
    t.append(["NIMH Controls","NA","NA","NA","NA","NA","v1.0"])                 
    t.append(["MGH IBD Prism Controls","NA","NA","NA","NA","NA","v1.1"])
    t.append(["iSAEC Controls","NA","NA","NA","NA","NA","v1.1"])
    t.append(["1000 Genomes","NA","NA","NA","NA","NA","v1.1 + Custom Content"])    
    makeTable(t,outFile)
    outFile.write("\n")
    outFile.write("\n")

    outFile.write( "Table 7: Merge Summary\n")
    t = [["nSamples", "nMale", "nFemale","nUnknownSex", "nEuropean","nSNPs (A+PA)", "nSNPs (X)", "nSNPs (Y)", "nSNPs (MT)", "nSNPs Total"]]
    root = inputVar["Output"] + inputVar["Project"] + "_CONTROLS"
    male, female, unknown = countMaleFemale(root + ".fam")
    apa, x, y, mt = countChr(root + ".bim")
    t.append([wc(root + ".fam"), male, female, unknown, uniq([inputVar["Output"] + "european.txt"]), apa, x, y, mt, wc(root + ".bim")])
    makeTable(t,outFile)
    outFile.write("\n")
    outFile.write("\n")

    outFile.write( "--------------------------------\n")

    outFile.write("\n")
    outFile.write("\n")
    outFile.write( "For Figure 1, see " + inputVar["Output"] + "euro.evec.ind.plot.png and " + inputVar["Output"] + "euro.evec.plot.pc12.png\n")


    return 0
