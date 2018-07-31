#! /usr/bin/python

import sys
import os

def writeMetaQC(inputVar, batches, batch):

    batchDir = inputVar["Output"] + batch + "/"
    if os.path.isdir(batchDir) == False:
        os.system("mkdir " + batchDir)

    metaFile = open(batchDir + "meta", 'w')

    metaFile.write("Script Directory: " + inputVar["Script"] + "\n")
    metaFile.write("Output Directory: " + inputVar["Output"] + batch + "/" + "\n")
    metaFile.write("Phenotype File: " + inputVar["Phenotype"] + "\n")        
    metaFile.write("Cohort File: " + inputVar["Cohort"] + "\n")
    metaFile.write("GAP Summary File: " + inputVar["GAP"] + "\n")
    metaFile.write("Project Name: " + batch + "\n")
    
    metaFile.write("gencall file: " + batches[batch]["qc"] + "\n")
    metaFile.write("zcall file: " + batches[batch]["qc"] + "\n")

    metaFile.write("gencall_missrate: " + batches[batch]["gencall_missrate"] + "\n")

    submit = 0
    
    metaFile.write("SAMPLE_QC: 0" + "\n")

    metaFile.write("RELATEDNESS_CHECK: 0" + "\n")

    metaFile.write("SEX_CHECK: 0" + "\n")

    metaFile.write("MENDEL_CHECK: 0" + "\n")
    
    if os.path.exists(batchDir + "pcaround1.checksum") == False:
        metaFile.write("PCA_ROUND1: 1" + "\n")
        submit += 1
    else:
        metaFile.write("PCA_ROUND1: 0" + "\n")

    metaFile.write("SITE_QC: 0" + "\n")        

    metaFile.close()

    if submit != 0:
        o = open(batchDir + "bad.samples.txt",'w')
        o.close()
    
        scriptDir = inputVar["Script"]
        q = batches[batch]["queue"]
        if q == "hour":
            q = "hour -W 4:00"
        
        c = "bsub -P " + batch + " -q " + q + " -o " + batchDir + "qc.out \"sh " + scriptDir + "scripts/qc.sh " +  batchDir + "meta" + "\""
        os.system(c)    

    return submit
