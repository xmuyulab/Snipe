import os
import sys
import pandas as pd
from time import time
import json
import pysam
import numpy as np
class snipeReOptions:
      path_core = ""
      read1 = ""
      read2 = ""
      id_file = ""
      out_path = ""
      threads = 1
      dict_ta = ""
      dict_te = "" 

'''
objective:Index the reference sequence of each strain separately,
          The index can be used directly when the filter threshold is 1.
input:
      -c:The folder where the core part of the reference sequence is located
      -1:The first sequence of samples
      -2:The second Sequence of Samples
      -o:Output catalogue 
      -idReport:The directory where the id report-file is located
      -t:Number of threads(default:t=1)  
      -dictTarget:Dictionary containing database information 
      -dictTemplate:Dictionary of statistical editing distance information

'''
#Corrector
def correct(sample_size,nm_sum):
    p1=2.2e-5
    p0=5e-10
    lambdan = np.power((1-p0)/(1-p1),sample_size-nm_sum)
    qr = np.power(p0/p1,nm_sum)
    return 1. / (1. + lambdan*qr)

def rec(snipeReOptions):
################### index_core model #################
    start = time()
    path = snipeReOptions.path_core
    r1 = snipeReOptions.read1
    r2 = snipeReOptions.read2
    t = snipeReOptions.threads
    o_path = snipeReOptions.out_path
    report = snipeReOptions.id_file
    dict_target = snipeReOptions.dict_ta
    dict_template = snipeReOptions.dict_te
    all_core = os.listdir(path)
    #if not os.path.exists('%s/coreindex'%(path)):
       #os.makedirs('%s/coreindex'%(path))

  #### index_core if index exists then continue##########
    for core in all_core:
        if os.path.exists('%s/%s/only_%s_blastn.1.bt2'%(path,core,core)):
           continue
        cmd1 = 'bowtie2-build  %s/%s/only_%s_blastn.fasta %s/%s/only_%s_blastn'%(path,core,core,path,core,core)
        os.system(cmd1)

  ###### rec  ######
    all_core = os.listdir(path)
    #Sample algin to special area of species
    ref_name = []
    counts0 = []
    counts1 = []
    counts2 = []
    f =open('%s'%(dict_target),'r')
    c = f.read()
    dict_target = eval(c)
    f.close()
    dict0 = {}
    for key in dict_target:
        if dict_target[key] not in dict0:
           dict0[dict_target[key]] = 0
    dict1 = dict0
    dict2 = dict0
    dict_sum = {}
    for core in all_core:
        cmd2 = 'bowtie2 -p %d -x %s/%s/only_%s_blastn -1 %s -2 %s -S %s/%s.sam'%(t,path,core,core,r1,r2,o_path,r1.split('/')[-1].split('.')[0])
        os.system(cmd2)
        print (cmd2)
        #if sam_file.endswith('.sam'):
        bf = pysam.AlignmentFile('%s/%s.sam'%(o_path,r1.split('/')[-1].split('.')[0]),'r')
        #else:
           #continue
        sc = core.split('_')[0] + ' '+core.split('_')[1]
        for r in bf:
            if r.cigarstring!=None:
                if 'S' not in str(r.cigarstring) and 'H' not in str(r.cigarstring):  
                    if r.has_tag("NM"):
                        if r.get_tag("NM") == 0:
                            dict0[sc]+=1
                        if r.get_tag("NM") == 1:
                            dict1[sc]+=1
                        if r.get_tag("NM") == 2:
                            dict2[sc]+=1
        for key in dict0:
            dict_sum[key] = dict0[key]+dict1[key]+dict2[key]
        bf.close()
        #os.remove('%s/%s.sam'%(o_path,core))
    dict_c = {}
    n_p = int(os.popen('less %s | grep @ | wc -l'%(r1)).read())
    for key in dict_sum: 
        dict_c[key] = correct(n_p,dict_sum[key])
    f_c =open('%s'%(r1.split('/')[-1].split('.')[0]),'w')
    f_c.write(str(dict_c))
    f_c.close()
    os.system('mv -f %s %s'%(r1.split('/')[-1].split('.')[0],o_path))
    cmd_c =  'cp %s %s/%s_temp.tsv'%(report,o_path,r1.split('/')[-1].split('.')[0])
    os.system(cmd_c)
    cmd_s = 'sed -i 1d %s/%s_temp.tsv'%(o_path,r1.split('/')[-1].split('.')[0])
    os.system(cmd_s)
    data = pd.read_csv('%s/%s_temp.tsv'%(o_path,r1.split('/')[-1].split('.')[0]),sep='\t')
    f =open('%s'%(dict_template),'r')
    c = f.read()
    dict_template = eval(c)
    f.close()
    abu = []
    abu_re = []
    da = pd.DataFrame()
    lst = data['Genome']
    lst1 = data['Final Best Hit Read Numbers']
    lst_genomes = []
    lst_nm0 = []
    lst_nmsum = []
    species = ''
    cor = []
    for g in lst:
        if g in dict_template:
           lst_genomes.append(dict_template[g])
        else:
           lst_genomes.append('')
        #print (dict_template[g])
        if len(g.split('_')) > 2:
           species = g.split('_')[0]+' '+g.split('_')[1]
        else:
           if g in dict_target:
              species = dict_target[g]
           if species in dict_c:
              cor.append(dict_c[species])
              lst_nm0.append(dict0[species])
              lst_nmsum.append(dict_sum[species])
           else:
              cor.append(0)
              lst_nmsum.append(0)
    for i in range(len(lst)):
        abu.append(lst1[i]/n_p)
        abu_re.append(lst1[i]*cor[i]/n_p)
    da['Genomes'] = lst_genomes
    da['Accession ID'] = lst
    da['Before Rectify'] = abu
    da['After Rectify'] = abu_re
    da['Final Guess'] = data['Final Guess']
    da['Final Best Hit'] = data['Final Best Hit']
    da['Final Best Hit Read Numbers'] = data['Final Best Hit Read Numbers']
    da['NM:sum'] = lst_nmsum
    da['Correction Probability']=cor
    da.to_csv('%s/%s_abundance.tsv'%(o_path,r1.split('/')[-1].split('.')[0]))
    os.remove('%s/%s_temp.tsv'%(o_path,r1.split('/')[-1].split('.')[0]))
    end = time()
    time_consume = end-start
    print ('rec time consuming: %d s'%(time_consume))
