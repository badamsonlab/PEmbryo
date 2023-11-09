##import matplotlib.pyplot as plt
##import seaborn as sns
import sys
import numpy as np
import os 
import time 
import pandas as pd
import gzip
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
np.set_printoptions(threshold=np.inf)

def update_Ns(reads,window):
    updated_reads=[]
    phred_33_scoring = list('!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ')
    for read in reads:
        r=""
        n_count=0
        seq = read[1]
        base_quals = read[3]
        base_scores = [phred_33_scoring.index(b) for b in base_quals]
        for b,base in enumerate(base_scores):
            if base < 30:
                r+="N"
            else:
                r+=seq[b]
        read[1]=r.upper()
        
        N_percentage = round(read[1][window[0]:window[1]].count("N")/(window[1]-window[0])*100,1)
        
        if N_percentage>10:
            pass
        else:
            updated_reads.append(read)
        
    return(updated_reads)


def validate_reads(reads,WT_ref):
    valid_reads = []

    for read in reads:
        matches=0
        seq = read[1]
            
        for i in range(4,4+40):
            if seq[i].upper()==WT_ref[i].upper() or seq[i].upper()=="N":
                matches+=1
        if round(matches/40*100,1)>=90:
            valid_reads.append(read)
        else:
            pass

    return(valid_reads)

match_dic = {}
match_dic[('A', 'A')] = 1
match_dic[('A', 'C')] = 0
match_dic[('A', 'T')] = 0
match_dic[('A', 'G')] = 0
match_dic[('C', 'A')] = 0
match_dic[('C', 'C')] = 1
match_dic[('C', 'T')] = 0
match_dic[('C', 'G')] = 0
match_dic[('G', 'A')] = 0
match_dic[('G', 'C')] = 0
match_dic[('G', 'T')] = 0
match_dic[('G', 'G')] = 1
match_dic[('T', 'A')] = 0
match_dic[('T', 'C')] = 0
match_dic[('T', 'T')] = 1
match_dic[('T', 'G')] = 0
match_dic[('N', 'N')] = 1
match_dic[('N', 'A')] = 1
match_dic[('N', 'T')] = 1
match_dic[('N', 'G')] = 1
match_dic[('N', 'C')] = 1
match_dic[('A', 'N')] = 1
match_dic[('T', 'N')] = 1
match_dic[('G', 'N')] = 1
match_dic[('C', 'N')] = 1

match_dic[('K', 'G')] = 1
match_dic[('K', 'T')] = 1
match_dic[('K', 'A')] = 0
match_dic[('K', 'C')] = 0
match_dic[('K', 'N')] = 1
match_dic[('G', 'K')] = 1
match_dic[('T', 'K')] = 1
match_dic[('A', 'K')] = 0
match_dic[('C', 'K')] = 0
match_dic[('N', 'K')] = 1

#K = G or T

#alignment_score = pairwise2.align.globalds(query_seq.upper(), ref_seq.upper(), match_dic, -3, -1, one_alignment_only=True,score_only=True)

query = sys.argv[1] ### input argument is target;edit e.g. chd2;G>A
target = str(query).split(';')[0]
edit = str(query).split(';')[1]

print(query)

time_stamp =str(round(time.time(),1))

def write_log(list_of_lists):
    with open('PEmbryo_log_'+target+"_"+edit+"_"+time_stamp+'.txt','a') as x:
        for line in list_of_lists:
            x.write("\t".join([str(x) for x in line])+"\n")
    return()

with open('PEmbryo_log_'+target+"_"+edit+"_"+time_stamp+'.txt','w') as f:
    f.write('log initiated\n')

write_log([['index','filename','group','sample','stage','target','edit','method','method_abbrev','date','run','editor','pegRNA','dnMLH1','status','# reads','# filtered_reads','# valid_reads','percent_after_filter','percent_after_validate',\
        'len(valid_reads)','wt','precise_edit','unintended_edit','percent_wt','percent_precise','percent_error']])

## INPUTS ##
folder = 'input_libraries/'
sample_list = pd.read_excel('PEmbryo_sample_sheet_master.xlsx',sheet_name="samples") # the sample list
ref_df = pd.read_excel('PEmbryo_sample_sheet_master.xlsx',sheet_name="edits") # the reference sequences

## PROCESS ## 

for s,sample in enumerate(sample_list.iterrows()):
    ##if sample_range[0]<=s<sample_range[1]:
    
    if sample[1]['Target'].lower() == target.lower() and (sample[1]['Edit'].lower() == edit.lower() or sample[1]['Edit'].lower()=='none'):#and s>686:

        reads = []
        try:
            sample_name=sample[1]['Sample'].replace(" ","_")
        except AttributeError:
            sample_name='sample_x'
        status=""
        date = str(sample[1]['date'])
        run = str(sample[1]['Run'])
        editor = str(sample[1]['Editor'])
        pegRNA = str(sample[1]['pegRNA'])
        dnMLH1 = str(sample[1]['dnMLH1'])
        stage = str(sample[1]['Stage'])

        ID = 'PEmbryo'+"_"+sample[1]['i5']+"_"+sample[1]['i7']+"_"+sample[1]['custom']+"_"+str(sample[1]['date'])+"_R1.fastq.gz"
        print(ID)

        ##well = sample[1]['well']
        ##plate = sample[1]['plate']
        group = sample[1]['Description'].replace(" ","_")
        method = sample[1]['Method']
        method_abbrev = method[0:3] #3 letter identifier of editing strategy

        if method=='control':
            method_abbrev = 'ctrl'

        for file in os.listdir(folder):
            if ID == file:
                filename=file
                with gzip.open(folder+file,'rt') as f:
                    data = [x.rstrip() for x in f.readlines()]
                    reads = [data[i-3:i+1] for i,x in enumerate(data) if (i+1)%4==0]

                    
                ### ref info

                flag=0
                ref_row=None
                for e in ref_df.iterrows():
                    #print(e[1]['target'],target,e[1]['edit'],edit,method_abbrev,[x for x in e[1]['method'].split(',')])
                    if e[1]['target'].lower() == target.lower() and e[1]['edit'].lower() == edit.lower() and method_abbrev.lower() in [x.lower() for x in e[1]['method'].split(',')]:
                        ref_row = e
                        flag=1
                        break

                if flag==1:
                    #window = [int(x) for x in ref_row[1]['window'].split("-")]   
                    cut_position = int(ref_row[1]['cut position'])
                    strand = ref_row[1]['strand']
                    if strand=='+':
                        window = [cut_position-18,cut_position+22]
                    elif strand=='-':
                        window = [cut_position-22,cut_position+18]
                    
                    WT_ref = ref_row[1]['WT_ref']
                    PE_ref = ref_row[1]['PE_ref']
                    status='found'
                else:
                    status='not_found'
                    #print('editing strategy not found')
                    #do something
                    
        
        ### filter reads
        
                filtered_reads = update_Ns(reads,window)
                valid_reads = validate_reads(filtered_reads,WT_ref)

                try:
                    percent_after_filter = round(len(filtered_reads)/len(reads)*100,1)
                    percent_after_validate = round(len(valid_reads)/len(filtered_reads)*100,1)
                except ZeroDivisionError:
                    percent_after_filter=0
                    percent_after_validate=0
            
                #print(s,filename,sample_name,len(reads),len(filtered_reads),len(valid_reads),'%/%',percent_after_filter,percent_after_validate)
        
        ### align reads

                wt=0
                precise_edit=0
                unintended_edit=0
                wt_reads=[]
                pe_reads=[]
                unintended_edit_reads=[]

                if flag==1:

                    for r,read in enumerate(valid_reads):
                        alignment_score = pairwise2.align.globalds(read[1][window[0]:window[1]].upper(),WT_ref[window[0]:window[1]].upper(), match_dic, -3, -1, one_alignment_only=True,score_only=True)
                        if alignment_score==int(window[1]-window[0]):
                            wt+=1
                            wt_reads.append(read)
                        else:
                            alignment_score = pairwise2.align.globalds(read[1][window[0]:window[1]].upper(),PE_ref[window[0]:window[1]].upper(), match_dic, -3, -1, one_alignment_only=True,score_only=True)
                            if alignment_score==int(window[1]-window[0]):
                                precise_edit+=1
                                pe_reads.append(read)
                            else:
                                unintended_edit+=1
                                unintended_edit_reads.append(read)
        
        #if ID.split("_")[1] in ['GTTTCGGA']:
        #    print("~~CONTROL~~")
        
                try:
                    percent_wt = round(wt/(wt+precise_edit+unintended_edit)*100,1)
                    percent_precise = round(precise_edit/(wt+precise_edit+unintended_edit)*100,1)
                    percent_unintended = round(unintended_edit/(wt+precise_edit+unintended_edit)*100,1)
                except ZeroDivisionError:
                    percent_wt = 0
                    percent_precise = 0
                    percent_unintended = 0
        
                #print(len(valid_reads),wt,precise_edit,unintended_edit,'%/%/%',percent_wt,percent_precise,percent_unintended)

                with open('output_libraries/wt_reads/'+ID[:-9]+'_'+target+'_wt_reads.fastq','w') as x:
                    for read in wt_reads:
                        for line in read:
                            x.write(line+'\n')

                with open('output_libraries/pe_reads/'+ID[:-9]+'_'+target+'_pe_reads.fastq','w') as x:
                    for read in pe_reads:
                        for line in read:
                            x.write(line+'\n')

                with open('output_libraries/error_reads/'+ID[:-9]+'_'+target+'_error_reads.fastq','w') as x:
                    for read in unintended_edit_reads:
                        for line in read:
                            x.write(line+'\n')

                write_log([[s,filename,group,sample_name,stage,target,edit,method,method_abbrev,date,run,editor,pegRNA,dnMLH1,status,len(reads),len(filtered_reads),len(valid_reads),percent_after_filter,percent_after_validate,wt+precise_edit+unintended_edit,wt,precise_edit,unintended_edit,percent_wt,percent_precise,percent_unintended]])


                break
        #break

print('successfully completed')


    
