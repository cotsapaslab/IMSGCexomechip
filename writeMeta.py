#! /usr/bin/python

import sys
import os

def writeMeta(inputVar, batches, batch):

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
    metaFile.write("maxJobs: " + inputVar["maxJobs"] + "\n")
    
    metaFile.write("gencall file: " + batches[batch]["gencall"] + "\n")
    metaFile.write("zcall file: " + batches[batch]["zcall"] + "\n")
    metaFile.write("minHet: " + batches[batch]["minHet"] + "\n")
    metaFile.write("maxHet: " + batches[batch]["maxHet"] + "\n")
    metaFile.write("maxNoCall: " + batches[batch]["maxNoCall"] + "\n")
    metaFile.write("maxMendelError: " + batches[batch]["maxMendelError"] + "\n")
    metaFile.write("gencall_missrate: " + batches[batch]["gencall_missrate"] + "\n")    
    metaFile.write("zcall_missrate: " + batches[batch]["zcall_missrate"] + "\n")
    metaFile.write("removesexerror: " + batches[batch]["removesexerror"] + "\n")        
    metaFile.write("removemendelerror: " + batches[batch]["removemendelerror"] + "\n")
    metaFile.write("unrelated_output: " + batches[batch]["unrelated_output"] + "\n")
    metaFile.write("use1KGdef: " + batches[batch]["use1KGdef"] + "\n")     

    submit = 0
    if os.path.exists(batchDir + "sampleqc.checksum") == False:
        metaFile.write("SAMPLE_QC: 1" + "\n")
        submit += 1
    else:
        metaFile.write("SAMPLE_QC: 0" + "\n")

    if os.path.exists(batchDir + "relatedness_check.checksum") == False:
        metaFile.write("RELATEDNESS_CHECK: 1" + "\n")
        submit += 1
    else:
        metaFile.write("RELATEDNESS_CHECK: 0" + "\n")

    if os.path.exists(batchDir + "sexcheck.checksum") == False:
        metaFile.write("SEX_CHECK: 1" + "\n")
        submit += 1
    else:
        metaFile.write("SEX_CHECK: 0" + "\n")

    if os.path.exists(batchDir + "mendelcheck.checksum") == False:
        metaFile.write("MENDEL_CHECK: 1" + "\n")
        submit += 1
    else:
        metaFile.write("MENDEL_CHECK: 0" + "\n")

    if os.path.exists(batchDir + "pcaround1.checksum") == False:
        metaFile.write("PCA_ROUND1: 1" + "\n")
        submit += 1
    else:
        metaFile.write("PCA_ROUND1: 0" + "\n")

    if os.path.exists(batchDir + "snpqc.checksum") == False:
        metaFile.write("SITE_QC: 1" + "\n")
        submit += 1
    else:
        metaFile.write("SITE_QC: 0" + "\n")        

    metaFile.close()

    if submit != 0:
        scriptDir = inputVar["Script"]
        q = batches[batch]["queue"]
        if q == "hour":
            q = "hour -W 4:00"

        rusage = batches[batch]["rusage"]
        c = "bsub -P " + batch + " -q " + q + " -R \"rusage[mem=" + rusage + "]\" -o " + batchDir + "qc.out \"sh " + scriptDir + "scripts/qc.sh " +  batchDir + "meta" + "\""
        os.system(c)    

    return submit
