#! /usr/bin/python

import sys
import os
import datetime

def makeLogFile(inputVar, batches):
    OUTPUT_DIR = inputVar["Output"]
    LOG_DIR = OUTPUT_DIR
#    LOG_DIR = OUTPUT_DIR + "logs/"
    if os.path.isdir(LOG_DIR) == False:
        os.system("mkdir " + LOG_DIR)

    n = 1
    LOG = ""    
    while LOG == "":
        if os.path.exists(LOG_DIR + "pipeline." + str(n) + ".log") == False:
            LOG = LOG_DIR + "pipeline." + str(n) + ".log"
            break
        else:
            n += 1

    print "Starting Exome Chip QC Pipeline..."
    print "See", LOG, "for log file."
    out = open(LOG, 'w')
    out.write("----------------------------" + "\n")
    out.write("Exome Chip QC Pipeline v4.0" + "\n")
    out.write("" + "\n")
    out.write(" ".join(["Meta File Used:", inputVar["Meta"]]) + "\n")
    out.write(" ".join(["Output Directory Used:", inputVar["Output"]]) + "\n")
    out.write(" ".join(["Script Directory Used:", inputVar["Script"]]) + "\n")
    out.write("" + "\n")

    for b in batches:
        q = batches[b]
        out.write("Parameters being used for " + b + ":" + "\n")
        for i in q:
            out.write(": ".join([i, q[i]]) + "\n")
        out.write(" " + "\n")

    out.write("----------------------------" + "\n")
    out.write(" " + "\n")    
    time = datetime.datetime.now()
    out.write("Analysis started at: " + time.strftime("%B %d %Y %H:%M:%S") + "\n")

    return out
    
        
