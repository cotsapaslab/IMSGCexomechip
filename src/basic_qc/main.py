#! /usr/bin/python

import sys
import os
import time

META = sys.argv[1] # input text file with all of input info

# Parse Meta File and initialize batch dictionary
qc_only = 0
skip_hla = 0
batches = {}
inputVar = {"Output":"","Phenotype":"","Cohort":"",
            "GAP":"","Project":"","Script":"","Meta":META,"maxJobs":"200","mergeMem":"4","mergeQueue":"hour -W 4:00"}
for line in open(META, 'r'):
    fields = line.split()
    if line.find("Output Directory:") != -1:
        OUTPUT_DIR = fields[2]
        if os.path.exists(OUTPUT_DIR) == False:
            print "Output directory does not exist. Exiting pipeline."
            sys.exit()
        if OUTPUT_DIR[len(OUTPUT_DIR)-1] != "/":
            OUTPUT_DIR += "/"
        inputVar["Output"] = OUTPUT_DIR

    elif line.find("Phenotype File:") != -1:
        PHENO = fields[2]
        if os.path.exists(PHENO) == False:
            print "Phenotype file does not exist. Exiting pipeline."
            sys.exit()        
        inputVar["Phenotype"] = PHENO
        
    elif line.find("Cohort File:") != -1:
        COHORT = fields[2]
        if os.path.exists(COHORT) == False:
            print "Cohort file does not exist. Exiting pipeline."
            sys.exit()        
        inputVar["Cohort"] = COHORT

    elif line.find("GAP Summary File:") != -1:
        GAP_SUMMARY = fields[3]
        inputVar["GAP"] = GAP_SUMMARY
        if os.path.exists(GAP_SUMMARY) == False:
            print "GAP summary file does not exist. Exiting pipeline."
            sys.exit()        

    elif line.find("Project Name:") != -1:
        PROJECT_NAME = fields[2:]
        if len(PROJECT_NAME) != 1:
            print "Incompatible project name: make sure it's one word with no spaces. Exiting pipeline."
            sys.exit()
        inputVar["Project"] = PROJECT_NAME[0]
        PROJECT_NAME = "".join(PROJECT_NAME)

    elif line.find("Script Directory:") != -1:        
        wd = fields[2]
        if os.path.exists(wd) == False:
            print "Script directory does not exist. Exiting pipeline."
            sys.exit()
        if wd[len(wd)-1] != "/":
            wd += "/"            
        inputVar["Script"] = wd

    elif line.find("QC_ONLY") != -1:
        qc_only = int(fields[1])
        inputVar["qc_only"] = qc_only

    elif line.find("SKIP_HLA") != -1:
        skip_hla = int(fields[1])
        inputVar["skip_hla"] = skip_hla

    elif line.find("maxJobsAtOneTime") != -1:
        inputVar["maxJobs"] = fields[1]

    elif line.find("mergeQueue") != -1:
        inputVar["mergeQueue"] = fields[1]

    elif line.find("mergeMem") != -1:
        inputVar["mergeMem"] = fields[1]

    elif len(fields) == 3:
        if fields[1] == "gencall" or fields[1] == "zcall" or fields[1] == "qc":
            batch = fields[0]            
            if batch not in batches:
                batches[batch] = {}

if inputVar["Script"] == "":
    print "Script Directory is not defined in META file. Exiting Pipeline."
    sys.exit()
if inputVar["Output"] == "":
    print "Output Directory is not defined in META file. Exiting pipeline."
    sys.exit()
if inputVar["Project"] == "":
    print "Project Name is not defined in META file. Exiting pipeline"
    sys.exit()
if len(batches) == 0:
    print "Define at least one input file in META file. Exiting pipeline"
    sys.exit()
    
# Make sure script directories exist
scrd = wd + "scripts/"
eigd = wd + "eigenstrat/"
hlad = wd + "hlaImputation/"
cond = wd + "controls/"
supd = wd + "suppl/"
D = [scrd,eigd,hlad,cond,supd]
exit = 0
for d in D:
    if os.path.exists(d) == False:
        print d, "does not exist. Exiting pipeline."
        exit = 1
if exit == 1:
    sys.exit()

# Add script directory to search path
sys.path.append(scrd)

# Import functions
from makeLogFile import *
from writeMeta import *
from writeMetaQC import *
from writeSuppMethods import *

# Parse input variables by batch
for line in open(META, 'r'):
    fields = line.split()
    if len(fields) == 3:
        batch = fields[0]
            
        if batch in batches:
            key = fields[1]
            value = fields[2]
            batches[batch][key] = value

# Fill in default parameters if not defined in meta file
Q = {"minHet":"0.1","maxHet":"1.20","gencall_missrate":"0.03",
     "zcall_missrate":"0.01","removesexerror":"1","removemendelerror":"1",
     "unrelated_output":"1","gencall":"","zcall":"","qc":"","skip":"0","queue":"hour","maxNoCall":"0.75","use1KGdef":"1",
     "maxMendelError":"10","rusage":"2"}


for batch in batches:
    for q in Q:
        if q not in batches[batch]:
            batches[batch][q] = Q[q]
        
for batch in batches:
    if batches[batch]["qc"] == "" and batches[batch]["gencall"] == "" and batches[batch]["zcall"] == "":
        print "Error in", batch, ". No input files defined"
        sys.exit()
    
    elif batches[batch]["qc"] != "" and batches[batch]["gencall"] != "" and batches[batch]["zcall"] != "":
        print "Error in", batch, ". Can only define either a set of gencall/zcall files or one qc files per batch. Exiting pipeline"
        sys.exit()

    elif batches[batch]["gencall"] != "" and batches[batch]["zcall"] == "":
        print "Error in", batch, ". Must define a zCall input."
        sys.exit()

    elif batches[batch]["gencall"] == "" and batches[batch]["zcall"] != "":
        print "Error in", batch, ". Must define a GenCall input."
        sys.exit()

    else:
        q = ["gencall","zcall","qc"]
        for i in q:
            f = batches[batch][i]
            if f != "":
                if os.path.exists(f + ".bed") == False:
                    print "Input file", f + ".bed", "for batch", batch,"does not exist. Exiting pipeline."
                    sys.exit()
                
# Make Log file
LOG = makeLogFile(inputVar, batches)

# Run basic QC on each batch separately than wait for processes to complete. If file QC'd already, only do prelim PCA to identify european samples
checkLOGS = []
for batch in batches:
    if batches[batch]["qc"] == "" and batches[batch]["skip"] == "0":
        bsubOut = OUTPUT_DIR + batch + "/qc.out"
        if os.path.exists(bsubOut) == True:
            os.system("rm " + bsubOut)
        p = writeMeta(inputVar,batches, batch)
        if p != 0:
            checkLOGS.append(bsubOut)            
        LOG.write(" ".join(["QC'ing", batch]) + "\n")
    elif batches[batch]["qc"] != "":
        bsubOut = OUTPUT_DIR + batch + "/qc.out"
        if os.path.exists(bsubOut) == True:
            os.system("rm " + bsubOut)
        p = writeMetaQC(inputVar,batches,batch)
        if p != 0:
            checkLOGS.append(bsubOut)                   
            LOG.write(" ".join(["Defining European Samples", batch]) + "\n")
    else:
        LOG.write(" ".join(["Skipping the QC of", batch]) + "\n")        


# Check outputs for each batch before merging
allSuccess = 1
while len(checkLOGS) != 0:
    for c in checkLOGS:
        if os.path.exists(c) == True:
            checkLOGS.pop(checkLOGS.index(c))
            success = 0
            for line in open(c, 'r'):
                if line.find("Successfully") != -1:
                    success = 1
            if success == 0:
                allSuccess = 0
                LOG.write(" ".join(["Error in bsub log file", c, "\n"]))            
    time.sleep(300)

if allSuccess == 0:
    LOG.write("Error in at least one batch... Exiting pipeline.\n")
    sys.exit()

exit = 0
for batch in batches:
    if batches[batch]["qc"] == "" and batches[batch]["skip"] == "0":
        final = OUTPUT_DIR + batch + "/" + batch
        if os.path.exists(final+".bed") == False:
            LOG.write("Error in " + batch + "\n")
            exit = 1

if exit == 1:
    LOG.write("Final batch files not present for at least one batch... Exiting pipeline.\n")    
    sys.exit()

if qc_only == 1:
    LOG.write("QC_ONLY option selected. Done with pipeline.\n")
    sys.exit()

# Make European samples file
c = "cat " + wd + "controls/european.txt > " + OUTPUT_DIR + "european.txt"
os.system(c)
for batch in batches:
    if batches[batch]["skip"] == "0":
        c = "cat " + OUTPUT_DIR + batch + "/" + "european.txt >> " + OUTPUT_DIR + "european.txt"
        os.system(c)

# Merge all batches with controls
if os.path.exists(OUTPUT_DIR + PROJECT_NAME + "_CONTROLS.bed") == False:
    c = "bsub -P pipeline -o merge.out -q " + inputVar["mergeQueue"] + " -R \"rusage[mem=" + inputVar["mergeMem"] + "]\" \"sh " + scrd + "mergeBatchesWControls.sh"
    q = [c]
    q.append(OUTPUT_DIR)
    q.append(PROJECT_NAME)
    q.append(wd)
    for batch in batches:
        if batches[batch]["qc"] != "" and batches[batch]["skip"] == "0":
            q.append(batches[batch]["qc"])
        elif batches[batch]["gencall"] != "" and batches[batch]["skip"] == "0":
            q.append(OUTPUT_DIR + batch + "/" + batch)
    q.append("\"")
    os.system(" ".join(q))
    checkLOGS = ["merge.out"]
    LOG.write(" ".join(["Merged batches together with controls", batch]) + "\n")


# Check outputs for each batch before merging
allSuccess = 1
while len(checkLOGS) != 0:
    for c in checkLOGS:
        if os.path.exists(c) == True:
            checkLOGS.pop(checkLOGS.index(c))
            success = 0
            for line in open(c, 'r'):
                if line.find("Successfully") != -1:
                    success = 1
            if success == 0:
                allSuccess = 0
    time.sleep(300)

if allSuccess == 0:
    LOG.write("Error in merging... Exiting pipeline.\n")
    sys.exit()


# Make merged cohorts file
c = "cat " + wd + "controls/cohorts.txt > " + OUTPUT_DIR + "cohorts.txt"
os.system(c)
for batch in batches:
    if batches[batch]["skip"] == "0":
        c = "cat " + OUTPUT_DIR + batch + "/" + "cohort.txt >> " + OUTPUT_DIR + "cohorts.txt"
        os.system(c)

if os.path.isdir(OUTPUT_DIR + "bsubLogs/") == True:
    if len(os.listdir(OUTPUT_DIR + "bsubLogs/")) != 0:
        c = "rm " + OUTPUT_DIR + "bsubLogs/*"
        os.system(c)

if os.path.exists(OUTPUT_DIR + "pca2.out") == True:
    os.system("rm " + OUTPUT_DIR + "pca2.out")
if os.path.exists(OUTPUT_DIR + "hla.out") == True:
    os.system("rm " + OUTPUT_DIR + "hla.out")
if os.path.exists(OUTPUT_DIR + "relatedness_check2.out") == True:
    os.system("rm " + OUTPUT_DIR + "relatedness_check2.out")
    
checkLOGS = []

# PCA
if os.path.exists(OUTPUT_DIR + "pcaround2.checksum") == False:
    q = "bsub -P " + PROJECT_NAME + " -q hour -W 4:00 -o " + OUTPUT_DIR + "pca2.out \"sh " + scrd + "pcaRound2.sh " + wd + " " + OUTPUT_DIR + " " + PROJECT_NAME + " " + OUTPUT_DIR + "cohorts.txt" + "\""
    os.system(q)
    checkLOGS.append(OUTPUT_DIR + "pca2.out")

# Check for related samples and dups after merging batches together
if os.path.exists(OUTPUT_DIR + "relatedness_check2.checksum") == False:
    q = "bsub -P " + PROJECT_NAME + " -q hour -W 4:00 -o " + OUTPUT_DIR + "relatedness_check2.out \"sh " + scrd + "relatednessCheck2.sh " + wd + " " + OUTPUT_DIR + " " + PROJECT_NAME  + " " + inputVar["maxJobs"] + "\""
    os.system(q)
    checkLOGS.append(OUTPUT_DIR + "relatedness_check2.out")

# HLA Imputation
if os.path.exists(OUTPUT_DIR + "hla.checksum") == False and skip_hla == 0:
    q = "bsub -P " + PROJECT_NAME + " -q week -o " + OUTPUT_DIR + "hla.out \"sh " + scrd + "hlaImputation.sh " + wd + " " + OUTPUT_DIR + " " + PROJECT_NAME + "\""
    os.system(q)
    checkLOGS.append(OUTPUT_DIR + "hla.out")


# Check bsub outputs before continuing
allSuccess = 1
while len(checkLOGS) != 0:
    for c in checkLOGS:
        if os.path.exists(c) == True:
            checkLOGS.pop(checkLOGS.index(c))
            success = 0
            for line in open(c, 'r'):
                if line.find("Successfully") != -1:
                    success = 1
            if success == 0:
                allSuccess = 0
                LOG.write(" ".join(["Error in bsub log file", c, "\n"]))            
    time.sleep(300)


if allSuccess == 0:
    LOG.write("Error in either HLA Imputation or PCA or relatedness check 2 ... Exiting pipeline.\n")
    sys.exit()

if os.path.exists(OUTPUT_DIR + "covar.euro.txt") == False:
    LOG.write("Error in PCA. No covar.euro.txt file completed.\n")
    sys.exit()

if os.path.exists(OUTPUT_DIR + PROJECT_NAME + "_CONTROLS.HLA.dosage") == False and skip_hla == 0:
    LOG.write("Error in HLA imputation. No .HLA.dosage file completed.\n")
    sys.exit()

s = 1
if os.path.exists(OUTPUT_DIR + "suppMethods.txt") == False:
    s = writeSuppMethods(inputVar,batches)
    if s != 0:
        LOG.write("Error in Writing up Supp Methods. Exiting before cleanup step.\n")
        sys.exit()
    
# Cleanup
x = ["y*","zCall*", "zcall*", "tmp*","gencall*", "keep*", "par*", "overlap*", "merge.geno", "merge.snp", "merge.ind",
     "pheno.txt", "*.hh", "*.nosex", "convertf.log","bsubLogs/","piHatCount.txt","meta"]

y = ["tmp*", "par*", "overlap*", "merge.geno", "merge.snp", "merge.ind","convertf.log"]

if os.path.exists(OUTPUT_DIR + "cleanup.checksum") == False:
    for b in batches:
        batchDir = OUTPUT_DIR + b + "/"
        if batches[b]["qc"] != "" and batches[b]["skip"] == "0":
            for i in y:
                w = "rm -R " + batchDir + i
                os.system(w)            
        elif batches[b]["qc"] == "" and batches[b]["skip"] == "0":
            if os.path.isdir(batchDir + "logs/") == False:
                w = "mkdir " + batchDir + "logs/"
                os.system(w)
            w = "mv " + batchDir + "*.log " + batchDir + "logs/"
            os.system(w)
            for i in x:
                w = "rm -R " + batchDir + i
                os.system(w)
            if os.path.exists(batchDir + "genome.genome") == True:
                w = "gzip " + batchDir + "genome.genome"
                os.system(w)


    if os.path.isdir(OUTPUT_DIR + "bsubLogs/") == True:
        w = "rm -R " + OUTPUT_DIR + "bsubLogs/"
        os.system(w)

    if skip_hla == 0:
        if os.path.exists(OUTPUT_DIR + inputVar["Project"] + "_CONTROLS.tmp.1.HLA.bim") == True:
            w = "rm -R " + OUTPUT_DIR + "*.tmp.*.HLA.*"
            os.system(w)

    if qc_only == 0:
        if os.path.exists(OUTPUT_DIR + "pcKeep0.fam") == True:
            w = "rm " + OUTPUT_DIR + "pcKeep*"
            os.system(w)
        if os.path.exists(OUTPUT_DIR + "pckeep.fam") == True:
            w = "rm " + OUTPUT_DIR + "pckeep*"
            os.system(w)            
        if os.path.exists(OUTPUT_DIR + "tmp.freq.frq") == True:
            w = "rm " + OUTPUT_DIR + "tmp.freq*"
            os.system(w)            
        if os.path.exists(OUTPUT_DIR + "tmp.bed") == True:
            w = "rm " + OUTPUT_DIR + "tmp*"
            os.system(w)
        if os.path.exists(OUTPUT_DIR + "keepMergeSNPs.txt") == True:
            w = "rm " + OUTPUT_DIR + "keep*"
            os.system(w)
        if os.path.exists(OUTPUT_DIR + "mergeList.txt") == True:
            w = "rm " + OUTPUT_DIR + "mergeList.txt"
            os.system(w)
        if os.path.exists(OUTPUT_DIR + "par.smartpca") == True:
            w = "rm " + OUTPUT_DIR + "par*"
            os.system(w)
        if os.path.exists(OUTPUT_DIR + "MERGE.eur.bim") == True:
            w = "rm " + OUTPUT_DIR + "MERGE*"
            os.system(w)
        if os.path.exists(OUTPUT_DIR + "merge.out") == True:
            w = "rm " + OUTPUT_DIR + "merge.out"
            os.system(w)                        
        if os.path.exists(OUTPUT_DIR + "convertf.log") == True:
            w = "rm " + OUTPUT_DIR + "convertf.log"
            os.system(w)
        if os.path.exists(OUTPUT_DIR + "smartpca.log") == True:
            w = "rm " + OUTPUT_DIR + "smartpca.log"
            os.system(w)                        
        if os.path.exists(OUTPUT_DIR + "genome.genome") == True:
            w = "gzip " + OUTPUT_DIR + "genome.genome"
            os.system(w)
        w = "rm *.hh"
        os.system(w)
        w = "rm *.nosex"
        os.system(w)
        w = "rm tmp*"
        os.system(w)
        

    w = "echo \"success\" > " + OUTPUT_DIR + "cleanup.checksum"
    os.system(w)
    
LOG.write("Done with pipeline\n")



