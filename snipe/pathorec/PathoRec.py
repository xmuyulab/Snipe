import os,sys,json
import pandas as pd
from time import time
import pysam
import numpy as np
class PathoRecOptions:
      path_core = ""
      read1 = ""
      read2 = ""
      id_file = ""
      out_path = ""
      threads = 1
      l_value = 1.0
      p0_value = 5e-10
      p1_value = 2.2e-5
      exp_tag = ""
      dict_ta = ""
      dict_te = "" 

'''
objective:Index the reference sequence of each strain separately,
          The index can be used directly when the filter threshold is 1.
input:
      -ssrRef:The folder where the core part of the reference sequence is located
      -1:The first sequence of samples
      -2:The second Sequence of Samples
      -o:Output catalogue 
      -idReport:The directory where the id report-file is located
      -t:Number of threads(default:t=1)  
      -dictTarget:Dictionary containing database information 
      -dictTemplate:Dictionary of statistical editing distance information

'''
# Rectification
def correct(l_v,p0_v,p1_v,sample_size,nm_sum):
    l = l_v
    p0 = p0_v
    p1 = p1_v
    lambdan = np.power((1-p0)/(1-p1),sample_size-nm_sum)
    qr = np.power(p0/p1,nm_sum)
    return 1. / (1. + l*lambdan*qr)

def rec(PathoRecOptions):
################### index_core model #################
    start = time()
    path = PathoRecOptions.path_core
    r1 = PathoRecOptions.read1
    r2 = PathoRecOptions.read2
    t = PathoRecOptions.threads
    l = PathoRecOptions.l_value
    p0 = PathoRecOptions.p0_value
    p1 = PathoRecOptions.p1_value
    o_path = PathoRecOptions.out_path
    report = PathoRecOptions.id_file
    dict_target = PathoRecOptions.dict_ta
    dict_template = PathoRecOptions.dict_te
    tag_name = PathoRecOptions.exp_tag
    all_core = os.listdir(path)

  #### index_core if index exists then continue##########
    for core in all_core:
        if os.path.exists('%s/%s/only_%s_blastn.1.bt2'%(path,core,core)):
           continue
        cmd1 = 'bowtie2-build  %s/%s/onlt_%s_blastn.fasta %s/%s/only_%s_blastn'%(path,core,core,path,core,core)
        os.system(cmd1)

  ###### rec  ######
    all_core = os.listdir(path)
    #Sample align to species specific regions
    ref_name = []
    counts0 = []
    counts1 = []
    counts2 = []
    f = open('%s'%(dict_target),'r')
    c = f.read()
    dict_target = eval(c)
    f.close()
    dict0 = {}
    dict1 = {}
    dict2 = {}
    dict_sum = {} 
    for key in dict_target:
        if dict_target[key] not in dict0:
           dict0[dict_target[key]] = 0
           dict1[dict_target[key]] = 0
           dict2[dict_target[key]] = 0
    for core in all_core:
        cmd2 = 'bowtie2 -p %d -x %s/%s/only_%s_blastn -1 %s -2 %s -S %s/%s_%s.sam'%(t,path,core,core,r1,r2,o_path,core,tag_name)
        os.system(cmd2)
        bf = pysam.AlignmentFile('%s/%s_%s.sam'%(o_path,core,tag_name),'r')
        sc = core.split('_')[0] + ' '+core.split('_')[1]
        for r in bf:
            if r.cigarstring != None:
                if 'S' not in str(r.cigarstring) and 'H' not in str(r.cigarstring):  
                    if r.has_tag("NM") and r.query_length > 50:
                        if r.get_tag("NM") == 0:
                            dict0[sc]+=1
                        if r.get_tag("NM") == 1:
                            dict1[sc]+=1
                        if r.get_tag("NM") == 2:
                            dict2[sc]+=1
     
        bf.close()
        os.remove('%s/%s_%s.sam'%(o_path,core,tag_name))
    for key in dict0:
        dict_sum[key] = dict0[key]+dict1[key]+dict2[key]
    dict_c = {}
    n_p = int(os.popen('less %s | grep @ | wc -l'%(r1)).read())
    for key in dict_sum: 
        dict_c[key] = correct(l,p0,p1,n_p,dict_sum[key])
    cmd_c =  'cp %s %s/%s_temp.tsv'%(report,o_path,tag_name)
    os.system(cmd_c)
    rebk = os.system("less %s/%s_temp.tsv | grep 'Total Number of Aligned Reads'"%(o_path,tag_name))
    if rebk==0:
        cmd_s = 'sed -i 1d %s/%s_temp.tsv'%(o_path,tag_name)
        os.system(cmd_s)
    data = pd.read_csv('%s/%s_temp.tsv'%(o_path,tag_name),sep='\t')
    f =open('%s'%(dict_template),'r')
    c = f.read()
    dict_template = eval(c)
    f.close()
    abu = []
    abu_rec = []
    final_rec = []
    da = pd.DataFrame()
    lst = data['Genome']
    lst1 = data['Final Best Hit Read Numbers']
    lst2 = data['Final Guess']
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
        abu_rec.append(lst1[i]*cor[i]/n_p)
        final_rec.append(lst2[i]*cor[i])
    da['Genomes'] = lst_genomes
    da['Accession ID'] = lst
    da['Rectified Final Guess'] = final_rec
    da['Final Guess'] = data['Final Guess']
    da['Rectified Probability']=cor
    da['SSR Aligned Reads'] = lst_nmsum
    da['Rectified Abundance'] = abu_rec
    da['Initial Abundance'] = abu
    da['Final Best Hit'] = data['Final Best Hit']
    da['Final Best Hit Read Numbers'] = data['Final Best Hit Read Numbers']
    da = da.sort_values(by=['Rectified Final Guess'],ascending=False)
    da.to_csv('%s/%s-rectified-report.tsv'%(o_path,tag_name),sep='\t',index=None)
    os.remove('%s/%s_temp.tsv'%(o_path,tag_name))
    end = time()
    time_consume = end-start
    print ('REC time consuming: %d s'%(time_consume))