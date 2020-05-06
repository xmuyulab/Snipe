import os
import sys
from time import time
class snipeIdOptions:
    ali_file = ""
    outdir = ""
    tag_name = ""

def id(snipeIdOptions):
    ################### snipeid model #################
    start =time()
    ali_file = snipeIdOptions.ali_file
    outdir = snipeIdOptions.outdir
    tag_n =snipeIdOptions.tag_name

    cmd3 = 'python2 pathoscope/pathoscope2.py ID -alignFile %s -fileType sam -outDir %s -expTag %s'%(ali_file,outdir,tag_n)
    os.system(cmd3)
    end = time()
    time_consume = end-start
    print ('snipeid time consuming: %d s'%(time_consume))
