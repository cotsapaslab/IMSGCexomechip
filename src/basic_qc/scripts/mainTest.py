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
            "GAP":"","Project":"","Script":"","Meta":META,"qc_only":"0","skip_hla":"0"}
for line in open(META, 'r'):
    fields = line.split()
    if line.find("Output Directory:") != -1:
        OUTPUT_DIR = fields[2]
        if os.path.exists(OUTPUT_DIR) == False:
            print "Output directory does not exist. Exiting pipeline."
            sys.exit()
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
        inputVar["Script"] = wd
        
    elif len(fields) == 3:
        if fields[1] == "gencall" or fields[1] == "zcall" or fields[1] == "qc":
            batch = fields[0]            
            if batch not in batches:
                batches[batch] = {}
    elif line.find("QC_ONLY") != -1:
        qc_only = int(fields[1])
        inputVar["qc_only"] = fields[1]

    elif line.find("SKIP_HLA") != -1:
        skip_hla = int(fields[1])
        inputVar["skip_hla"] = fields[1]
        
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
     "maxMendelError":"10"}


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
                

writeSuppMethods(inputVar,batches)


