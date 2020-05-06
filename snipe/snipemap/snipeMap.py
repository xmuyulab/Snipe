import os
import sys
import pandas as pd
from time import time
import json
import pysam
import numpy as np
class snipeOptions:
    n = 1
    read1 = ""
    read2 = ""
    reference_path = ""
    filter_path = ""
    core_path = ""
    tag_n = ""
    out_path = ""
    outAlignFile = "outalign.sam"
    threads = 1
    index_path="."
'''
objective:The filter compares the samples to the core parts of various strains.
          and then screens one or more strains as a reference sequence.
input:
     -1:The first sequence of samples
     -2:The second Sequence of Samples
     -referenceFiles:The directory where the reference sequence is located
     -filterFiles:The directory where the filter sequence is located
     -coreRefFiles:Catalogue of Core Parts of Strains
     -outDir:Output catalogue 
     -t:Number of threads(default:t=1)     
'''

def map(snipeOptions):
    start =time()
    ################### snipeFilter model #################
    r1 = snipeOptions.read1
    r2 = snipeOptions.read2
    ref_path = snipeOptions.reference_path
    fil_path = snipeOptions.filter_path
    o_path = snipeOptions.out_path
    tag_name = snipeOptions.tag_n
    i_p = snipeOptions.index_path
    outAlignF = snipeOptions.outAlignFile
    t = snipeOptions.threads
    cmd2 = 'python2 pathoscope/pathoscope2.py MAP -1 %s -2 %s -targetRefFiles %s -filterRefFiles %s -outDir %s -outAlign %s -indexDir %s -numThreads %d -expTag %s'%(r1,r2,ref_path,fil_path,o_path,outAlignF,i_p,t,tag_name)
    os.system(cmd2)
    end = time()
    time_consume = end-start
    print ('snipeMap time consuming: %d s'%(time_consume))
