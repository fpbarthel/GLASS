# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 13:10:57 2019

@author: medgnt
"""

import matplotlib.pyplot as plt 
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Extract shared and private VAFs from paired samples.")
parser.add_argument('-c', '--genotypecomparison', dest='genotypecomparison', required=True, type=str, help='Extracted GLASS_genotype_comparison file. E.g. GLASS_genotype_comparison_extracted.tsv')
parser.add_argument('-g', '--genotypes', dest='genotypes', required=True, type=str, help='GLASS_genotypes file. E.g. GLASS_genotypes.csv')
parser.add_argument('-s', '--set', dest='set', required=True, type=str, help='Set of tumours to include. E.g. silver_set.csv')
parser.add_argument('-t', '--titan', dest='titan', required=True, type=str, help='TITAN parameter estimates. E.g. titanparams_synapse.tsv')
parser.add_argument('-o', '--outdir', dest='outdir', required=True, type=str, help='Output directory')


args = parser.parse_args()

cpns={}        
with open(args['genotypecomparison'],'r') as file:    
    for line in file:
        l=line.strip().split('\t')
        if len(l)==6:  
            if l[0] not in cpns:
                cpns[l[0]]={}
            cpns[l[0]][str(l[1])+str(l[2])+str(l[3])]=[l[4],l[5]]
            

muts={} 
pairs={}      
with open(args['set'],'r') as file:    
    file.readline()    
    for line in file:
        l=line.strip().split(',')
        muts[l[1][:12]]={}
        muts[l[1]][l[2]]={}
        muts[l[1]][l[3]]={}
        pairs[l[0]]=[l[2],l[3]]
        
with open(args['genotypes'],'r') as file:
    file.readline()    
    for line in file:
        l=line.strip().split(',')
        if l[0][:12] in muts:
            if l[0] in muts[l[0][:12]]:
                #create a dictionary of {case:{sample:{mutation:[ref_counts,alt_counts,read_depth,mutect_called?]}}}
                muts[l[0][:12]][l[0]][str(l[1])+str(l[2])+str(l[4])]=[l[5],l[6],l[7],l[8]]

       
samples={}
with open(args['titan'],'r') as file:
    file.readline()
    for line in file:
        l=line.strip().split('\t')
        if l[0]!='':
            samples[l[2][1:19]]=[l[5],l[7]]

for case in muts:
    for sam in muts[case]:
        ploidy=float(samples[sam[:18]][1].strip('"'))
        for mut in muts[case][sam]:          
            #validate if mutect called
            if muts[case][sam][mut][3]=='t':
                muts[case][sam][mut].append('t')
            else:
                muts[case][sam][mut].append('f')
            
            #calculate coverage
            cov=int(muts[case][sam][mut][0])+int(muts[case][sam][mut][1])
            muts[case][sam][mut].append(cov)
    
            #mark as shared or private
            status='unknown' 
            for pair in pairs:
                if sam==pairs[pair][0]: #if sample is primary
                    other=pairs[pair][1]
                    break
                elif sam==pairs[pair][1]: #if sample is recurrent
                    other=pairs[pair][0]
                    break
            if mut in muts[case][other]:
                if muts[case][other][mut][3]=='t':
                    status='shared'
                elif int(muts[case][other][mut][0])+int(muts[case][other][mut][1])>=30:
                    status='private'
            muts[case][sam][mut].append(status)        
            
            #mark as copy neutral or cna region
            logr=''
            if mut in cpns[pair]:
                if sam==pairs[pair][0]: #if sample is primary
                    if cpns[pair][mut][0]!='':
                        logr=float(cpns[pair][mut][0])   
                elif sam==pairs[pair][1]: #if sample is recurrent
                    if cpns[pair][mut][1]!='':
                        logr=float(cpns[pair][mut][1])
                        
            if logr=='': 
                muts[case][sam][mut].append('-')
            else:    
                if round(ploidy-logr)==0:
                    muts[case][sam][mut].append('neutral')
                else:
                    muts[case][sam][mut].append('cna')
                
            #mark as clonal or subclonal
            if muts[case][sam][mut][1]=='0':
                muts[case][sam][mut].append('0')
            elif float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1]))<=0.24:
                muts[case][sam][mut].append('subclonal')
            else:
                muts[case][sam][mut].append('clonal')
                                
with open(args['outdir']+'metadata.tsv','w+') as file:
    file.write('case\tsample\tsharedclonal\tsharedsubclonal\tprivateclonal\tprivatesubclonal\tallclonal\tallsubclonal\tmeandepth\tpurity\tploidy\tminpurity\tmaxpurity\troundedploidy\n')
    for case in muts:
        for sam in muts[case]:
            sharedclonal=0
            sharedsubclonal=0
            privateclonal=0
            privatesubclonal=0
            allclonal=0
            allsubclonal=0
            sharedvafs=[]
            privatevafs=[]
            allvafs=[]
            count=0
            depths=[]
            purity=float(samples[sam[:18]][0].strip('"'))
            ploidy=float(samples[sam[:18]][1].strip('"'))
            for mut in muts[case][sam]:
                depths.append(int(muts[case][sam][mut][2]))
                count+=1
                if muts[case][sam][mut][4]=='t':
                    if int(muts[case][sam][mut][5])>=30:
                        if muts[case][sam][mut][7]=='neutral':
                            if muts[case][sam][mut][6]=='shared':
                                if muts[case][sam][mut][8]=='clonal':
                                    sharedclonal+=1
                                    sharedvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                                    allclonal+=1
                                    allvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                                elif muts[case][sam][mut][8]=='subclonal':
                                    sharedsubclonal+=1
                                    sharedvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                                    allsubclonal+=1
                                    allvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                            elif muts[case][sam][mut][6]=='private':
                                if muts[case][sam][mut][8]=='clonal':
                                    privateclonal+=1
                                    privatevafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                                    allclonal+=1
                                    allvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                                elif muts[case][sam][mut][8]=='subclonal':
                                    privatesubclonal+=1
                                    privatevafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                                    allsubclonal+=1
                                    allvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                            elif muts[case][sam][mut][6]=='unknown':
                                if muts[case][sam][mut][8]=='clonal':
                                    allclonal+=1
                                    allvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
                                elif muts[case][sam][mut][8]=='subclonal':
                                    allsubclonal+=1
                                    allvafs.append(float(muts[case][sam][mut][1])/(float(muts[case][sam][mut][0])+float(muts[case][sam][mut][1])))
            file.write(case+'\t'+sam+'\t'+str(sharedclonal)+'\t'+str(sharedsubclonal)+'\t'+str(privateclonal)+'\t'+str(privatesubclonal)+'\t'+str(allclonal)+'\t'+str(allsubclonal)+'\t'+str(sum(depths)/float(count))+'\t'+str(purity)+'\t'+str(ploidy)+'\t'+str(float(purity)-0.1)+'\t'+str(float(purity)+0.1)+'\t'+str(round(float(ploidy)))+'\n')
            
            with open(args['outdir']+sam+'_private.txt','w+') as file2:        
                for i in privatevafs:
                    file2.write(str(i)+'\n')
            plt.figure()
            plt.hist(privatevafs,bins=np.arange(0.0, 1.0, 0.01))
            plt.ylabel('Count')
            plt.xlabel('VAF')
            plt.xlim([0,1])
            plt.savefig(args['outdir']+sam+'_private.pdf')
            with open(args['outdir']+sam+'_shared.txt','w+') as file2:        
                for i in sharedvafs:
                    file2.write(str(i)+'\n')
            plt.figure()
            plt.hist(sharedvafs,bins=np.arange(0.0, 1.0, 0.01))
            plt.ylabel('Count')
            plt.xlabel('VAF')
            plt.xlim([0,1])
            plt.savefig(args['outdir']+sam+'_shared.pdf')
            with open(args['outdir']+sam+'_all.txt','w+') as file2:        
                for i in allvafs:
                    file2.write(str(i)+'\n')
            plt.figure()
            plt.hist(allvafs,bins=np.arange(0.0, 1.0, 0.01))
            plt.ylabel('Count')
            plt.xlabel('VAF')
            plt.xlim([0,1])
            plt.savefig(args['outdir']+sam+'_all.pdf')

