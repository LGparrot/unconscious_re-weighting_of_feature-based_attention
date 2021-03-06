#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:32:22 2018

@author: lasseguldener
"""
from __future__ import division
import pickle 
import numpy as np
import csv
import scipy 
import os
import math
import pandas as pd #makes it easy to store dictionaries in csv 
from scipy.stats import norm
from math import exp,sqrt
from astropy.stats import median_absolute_deviation
from scipy.stats import norm
from math import exp,sqrt
from datetime import datetime
### BEHAV anaylsis '''



'''change path here  '''

wdir=os.getcwd()
datapath=wdir + u'/Downloads/data/' 

datapath1=wdir+ u'/results/pre_traget/'


'''if sub is even: block is even = right-weighted, block = uneven left_weighted
if sub is uneven: block is even = left_weightes, block = uneven right_weighted'''  

nTrials=36

count_nan=0




nSubs=[1,2,3,4,7]#b=14
#nSubs=[8]#b=12
#nSubs=[9,10,11]#b=14
#nSubs=[14,15,16,17]#b=10
nSubs=[1,2,3,4,7,8,9,10,11,14,15,16,17]

#
nBlocks=14#[1,2,3,4,5,6,7,8,9,10,11,12,13,14]#14
co=3  #cut off for outlier removal in SD   
    
#this func finds the correct indices 
def all_indices(value, qlist):
        indices = []
        idx = -1
        while True:
            try:
                idx = qlist.index(value, idx+1)
                indices.append(idx)
            except ValueError:
                break
        return indices


''' this is how I try to avoid excessive loops to asign values to keys in dic according to indices'''

def TrialFinder_switch(ind,ind2,dic_entry,dv,block): #ind = indices (list), dic=dic with condis for odd and even , dv = dependent variable
    if block % 2==0:
        if len(ind)>0:
            av=np.zeros(shape=(len(ind)))
            av=av.tolist()        

            for i in range(len(ind)):
                av[i] = dv[ind[i]]
            
            dic_entry.append(av)
        
        else:
            av=[]
            dic_entry.append(av)

    elif block % 2!=0:
        if len(ind2)>0:
            av=np.zeros(shape=(len(ind2)))
            av=av.tolist()        

            for i in range(len(ind2)):
                av[i] = dv[ind2[i]]
            
            dic_entry.append(av)
        
        else:
            av=[]
            dic_entry.append(av)    
    
    return dic_entry


def DIC_WRITER(dic_list,labels): #put list with title and lables of the dictionary you want to create 
   for title in dic_list:
        title={}    
        for l in labels:#loop through lable, each key becomes an empty list
            title['{0}'.format(l)]=[]
        
        return title

Z = norm.ppf  
 
def dPrime(hits, misses, fas, crs):
    vals=[hits,misses,fas,crs]
    if all(vals)>0:
    # Floors an ceilings are replaced by half hits and half FA's
        halfHit = 0.5/(hits+misses)
        halfFa = 0.5/(fas+crs)
     
        # Calculate hitrate and avoid d' infinity
    
        hitRate = hits/(hits+misses)
        if hitRate == 1: hitRate = 1-halfHit
        if hitRate == 0: hitRate = halfHit
     
        # Calculate false alarm rate and avoid d' infinity
    
        faRate = fas/(fas+crs)
        if faRate == 1: faRate = 1-halfFa
        if faRate == 0: faRate = halfFa
     
        # Return d', beta, c and Ad'
        out = {}
        out['d'] = Z(hitRate) - Z(faRate) #that´s for a yes/no task!!!
        out['beta'] = exp((Z(faRate)**2 - Z(hitRate)**2)/2)
        out['c'] = -(Z(hitRate) + Z(faRate))/2
        out['Ad'] = norm.cdf(out['d']/sqrt(2))
        
    else:    
        out = {}
        out['d']=999
        out['beta']=999
        out['c'] =999 
        out['Ad'] =999
  
    return out   

#
def Beta(H,F):# this is for non parametric B'' 
    if H >0:
        if F >0:
            if H>=F:
                b=(H*(1-H)-F*(1-F))/(H*(1-H)+F*(1-F))
            elif H<F:
                b=(F*(1-F)-H*(1-H))/(F*(1-F)+ H*(1-H))
        else:b='nan'
    else: b='nan'        
    
    return b


def MAD(key):
    median=np.nanmedian(key)
    mini=median - 3*(median_absolute_deviation(key))
    maxi=median + 3*(median_absolute_deviation(key))
    key=[x for x in key if x > mini] 
    key=[x for x in key if x < maxi]
    


def trimm(list_n):
    if len(list_n)>0:
        mean=np.nanmean(list_n)  
        sd=np.nanstd(list_n)
        upper_limit=mean + co*sd
        lower_limit=mean - co*sd
        list_n=[x for x in list_n if x > lower_limit] 
        list_n=[x for x in list_n if x < upper_limit]
        list_n=np.nanmean(list_n)
    
    else:
        list_n.append('nan')
    
    return list_n

RT=[]
avRT=[]
SW=[]
AL=[]
SUB=[]
for sub in nSubs:
    if sub < 8:
        nBlocks=14
    elif sub ==8:
        nBlocks=12
    elif sub > 8 and sub < 12:
        nBlocks=14
    elif sub > 13:
        nBlocks=10


    #nBlocks=10
    
    '''make some dic for SA'''
    
    lables=['Hit1', 'Missed_Hit1','False_Alarm1', 'Correct_Reject1','Hit2', 'Missed_Hit2','False_Alarm2', 'Correct_Reject2',
            'Hit3', 'Missed_Hit3','False_Alarm3', 'Correct_Reject3','Hit4', 'Missed_Hit4','False_Alarm4', 'Correct_Reject4'] 
    
    SA_even={}
    for l in lables:
        SA_even['{0}'.format(l)]=[]
    
    SA_uneven={}
    for l in lables:
        SA_uneven['{0}'.format(l)]=[]

    ''' make dictionaries containing lists to collect the values of each condition across all blocks '''
    freq2infreq_rt={}
    infreq2freq_rt={}
    freq2vertic_rt={}
    infreq2vertic_rt={}
    vertic2freq_rt={}
    vertic2infreq_rt={}
    infreq2vertic_rt={}
    freq2vertic_rt={}
    ns_freq_rt={}
    ns_infreq_rt={}
    ns_vertic_rt={}
       
    ''' fill the dictionary, every entry becomes an empty list for data collection. 
    These can be flattened and used for mean calculation. '''
    
    lvls=[1,2,3,4]
    for l in lvls:
        freq2infreq_rt['{0}'.format(l)]=[]
        infreq2freq_rt['{0}'.format(l)]=[]
        freq2vertic_rt['{0}'.format(l)]=[]
        infreq2vertic_rt['{0}'.format(l)]=[]
        vertic2freq_rt['{0}'.format(l)]=[]
        vertic2infreq_rt['{0}'.format(l)]=[]
        ns_freq_rt['{0}'.format(l)]=[]
        ns_infreq_rt['{0}'.format(l)]=[]
        ns_vertic_rt['{0}'.format(l)]=[]
        infreq2vertic_rt['{0}'.format(l)]=[]
        freq2vertic_rt['{0}'.format(l)]=[]
        
    ''' dic to collect no. of trials for each condition'''
    
    '''switch'''
    n_freq2infreq={}
    n_infreq2freq={}
    n_freq2vertic={}
    n_infreq2vertic={}
    n_vertic2infreq={}
    n_vertic2freq={}
    '''no switch'''
    n_ns_freq={}
    n_ns_infreq={}
    n_ns_vertic={}
    
    for l in range(1,5):
        n_freq2infreq['{0}'.format(l)]=[] 
        n_infreq2freq['{0}'.format(l)]=[]     
        n_freq2vertic['{0}'.format(l)]=[] 
        n_infreq2vertic['{0}'.format(l)]=[]     
        n_vertic2freq['{0}'.format(l)]=[] 
        n_vertic2infreq['{0}'.format(l)]=[]  
        n_ns_vertic['{0}'.format(l)]=[]  
        n_ns_freq['{0}'.format(l)]=[]  
        n_ns_infreq['{0}'.format(l)]=[]  
    
    '''dictionray that will later contain number of trials for each level of awareness '''
  
    
    aware_lvls={}
    for a in range(1,5):
        aware_lvls['all_AL{0}'.format(a)]=[]
        
    all_nans=np.zeros(shape=(nBlocks))
    all_nans_recog=np.zeros(shape=(nBlocks))
    ''' some list to collect rates for sensitivity analysis '''
    H=[]
    M=[]
    F=[]
    C=[]
    
    '''make lists to collect values for regression format'''
    rt=[]
    oris=[]
    ALs=[]
    weighting=[]
    ID=[]
    o=[]
    all_blocks=[]
    
    all_actual={}
    all_preds={}
    
    for a in range(1,5):
        all_actual['{0}'.format(a)]=[]
        all_preds['{0}'.format(a)]=[]
        
    bcount=0#
    nancount=0
    for block in range(nBlocks):
        bcount=bcount+1
        awareness_levels = np.zeros(shape=(nTrials)) #will be overwritten by a consequent block
        sub_acc=np.zeros(shape=(nTrials))
        sub_rt=np.zeros(shape=(nTrials))
        ori=np.zeros(shape=(nTrials))
        dummy=np.zeros(shape=(nTrials))
        switch=np.zeros(shape=(nTrials))
        prior_ori=np.zeros(shape=(nTrials))

        '''lists for ROC analysis'''
        actual=np.zeros(shape=(nTrials))
        predictions=np.zeros(shape=(nTrials))
        
        file = open((datapath + 'mainexp_pilot_sub%.3i_block%.3i.pkl' %(sub,block)), 'rb')
        out_wmt = pickle.load(file,encoding='latin1')# this does the trick if you get trouble loading the pickel (also 'bytes' should work)
        file.close()
        all_blocks.append(out_wmt) # now all files in one list 
        
        count_nan=0 
        count_nan_recog=0
        
        for t in range(nTrials):#loop through file list
            ID.append(sub)
            if sub < 8:
                if sub % 2 ==0:
                    if block %2 ==0:
                        weighting.append('right')
                    else:weighting.append('left')
                    
                elif sub % 2 !=0:
                    if block %2 ==0:
                        weighting.append('left')
                    else:weighting.append('right')
            else:
                if sub % 2 !=0:
                    if block %2 ==0:
                        weighting.append('right')
                    else:weighting.append('left')
                    
                elif sub % 2 ==0:
                    if block %2 ==0:
                        weighting.append('left')
                    else:weighting.append('right')
                
            rt.append(all_blocks[block][t]['sub_rt'])   

            dummy[t]=1  
            sub_acc[t]=all_blocks[block][t]['sub_acc']
            

            if math.isnan(all_blocks[block][t]['sub_rt']):
                sub_rt[t]= 0

            else:
                sub_rt[t]=all_blocks[block][t]['sub_rt'] + 0.382


            
            if all_blocks[block][t]['sub_recog'] == 1: # AL1 awareness  scale point 1 = unaware
                awareness_levels[t]= 1
                ALs.append(1)
            elif all_blocks[block][t]['sub_recog'] == 2: # AL2 awareness  scale point 2 = maybe saw something
                awareness_levels[t]= 2  
                ALs.append(2)
            elif all_blocks[block][t]['sub_recog'] == 3: # AL3 awareness  scale point 3 = almost aware
                awareness_levels[t]= 3
                ALs.append(3)
            elif all_blocks[block][t]['sub_recog'] == 4: # AL4 awareness  scale point 4 = fully aware
                awareness_levels[t]= 4
                ALs.append(4)
            elif math.isnan(all_blocks[block][t]['sub_recog']):# if nan , any trial with missing value for aware
                awareness_levels[t]= 0
                ALs.append('nan')
                count_nan_recog=count_nan_recog+1
            

                    

            if math.isnan(all_blocks[block][t]['sub_resp']):#count NaNs
                count_nan=count_nan + 1
                
            if all_blocks[block][t]['gabor_ori']== 180.0:               #getACCs and RTs for all Oris
                ori[t]=1 #vertic
                oris.append(1)
                predictions[t]=1
                if sub_acc[t]==1: #correctly identified as vertical (correctly rejecting it beeing non-vertical )
                    C.append(1)
                elif sub_acc[t]==1:
                    F.append(1)
           
            elif all_blocks[block][t]['gabor_ori']< 180.0:
                ori[t]=2 #left
                oris.append(2)
                predictions[t]=2
                if sub_acc[t]==1:
                    H.append(1)
                elif sub_acc[t]==0:
                    M.append(1)
                
            elif all_blocks[block][t]['gabor_ori'] > 180.0:
                ori[t]=3 #right
                oris.append(3)
                predictions[t]=2
                if sub_acc[t]==1:
                    H.append(1)
                elif sub_acc[t]==0:
                    M.append(1) # not recognized 
                    
            #if all_blocks[block][t-1]['gabor_ori']==all_blocks[block][t]['gabor_ori']:
            if ori[t]==ori[t-1]:    
                switch[t]=0
            else:switch[t]=1
            
            #here omit any pre-switch trials with AL2-4
            if t >0:
                if switch[t]==1:
                    if ALs[t-1]!=1:
                        sub_rt[t]=0 
                
            if all_blocks[block][t-1]['gabor_ori']==180:prior_ori[t]=1 #vertical
            elif all_blocks[block][t-1]['gabor_ori']<180:prior_ori[t]=2#left
            elif all_blocks[block][t-1]['gabor_ori']>180:prior_ori[t]=3#right
            actual[t]=all_blocks[block][t]['sub_resp']
            
            if sub_rt[t]==0:
                nancount=nancount+1
                sub_rt[t]='nan'

            all_nans_recog[block]=count_nan_recog
            all_nans[block]=count_nan
        
            if awareness_levels[t]==1:
                all_preds['1'].append(predictions[t])
                all_actual['1'].append(actual[t])
            elif awareness_levels[t]==2:
                all_preds['2'].append(predictions[t])
                all_actual['2'].append(actual[t])
            elif awareness_levels[t]==3:
                all_preds['3'].append(predictions[t])
                all_actual['3'].append(actual[t])
            elif awareness_levels[t]==4:
                all_preds['4'].append(predictions[t])
                all_actual['4'].append(actual[t])
                
        print (sub_rt)       
        condi_zip=list(zip(ori,awareness_levels))
        aware_zip=list(zip(dummy,awareness_levels))
        sa_zip=list(zip(ori,sub_acc,awareness_levels))
        condi_zip1=list(zip(ori,awareness_levels,switch))
        aware_zip1=list(zip(dummy,awareness_levels,switch))
        sa_zip1=list(zip(ori,sub_acc,awareness_levels,switch))
        
        switchzip=list(zip(ori,awareness_levels,switch,prior_ori))
        
        '''switch indices taking prior ori into account'''
        sw_right2vertic_AL1_ind=all_indices((1,1,1,3),switchzip)
        sw_right2vertic_AL2_ind=all_indices((1,2,1,3),switchzip)
        sw_right2vertic_AL3_ind=all_indices((1,3,1,3),switchzip)
        sw_right2vertic_AL4_ind=all_indices((1,4,1,3),switchzip)

        sw_right2left_AL1_ind=all_indices((2,1,1,3),switchzip)
        sw_right2left_AL2_ind=all_indices((2,2,1,3),switchzip)
        sw_right2left_AL3_ind=all_indices((2,3,1,3),switchzip)
        sw_right2left_AL4_ind=all_indices((2,4,1,3),switchzip)
        
        sw_left2vertic_AL1_ind=all_indices((1,1,1,2),switchzip)
        sw_left2vertic_AL2_ind=all_indices((1,2,1,2),switchzip)
        sw_left2vertic_AL3_ind=all_indices((1,3,1,2),switchzip)
        sw_left2vertic_AL4_ind=all_indices((1,4,1,2),switchzip)
                
        sw_left2right_AL1_ind=all_indices((3,1,1,2),switchzip)
        sw_left2right_AL2_ind=all_indices((3,2,1,2),switchzip)
        sw_left2right_AL3_ind=all_indices((3,3,1,2),switchzip)
        sw_left2right_AL4_ind=all_indices((3,4,1,2),switchzip)
        
        sw_vertic2right_AL1_ind=all_indices((3,1,1,1),switchzip)
        sw_vertic2right_AL2_ind=all_indices((3,2,1,1),switchzip)
        sw_vertic2right_AL3_ind=all_indices((3,3,1,1),switchzip)
        sw_vertic2right_AL4_ind=all_indices((3,4,1,1),switchzip) 
        
        sw_vertic2left_AL1_ind=all_indices((3,1,1,1),switchzip)
        sw_vertic2left_AL2_ind=all_indices((3,2,1,1),switchzip)
        sw_vertic2left_AL3_ind=all_indices((3,3,1,1),switchzip)
        sw_vertic2left_AL4_ind=all_indices((3,4,1,1),switchzip) 
        
        if sub <= 8:    #gerade blöcke rechtsgewichtet
            if sub  % 2 != 0:
                if bcount % 2 ==0:
                    
                    n_freq2infreq['1'].append(len(sw_right2left_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_right2left_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_right2left_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_left2right_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_left2right_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_left2right_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_left2right_AL4_ind))
                
                    n_freq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                   
                    n_infreq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
    
                    n_vertic2freq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2left_AL4_ind))
                    
                    
                    
    
                elif bcount % 2 !=0:
                    n_freq2infreq['1'].append(len(sw_left2right_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_left2right_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_left2right_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_left2right_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_right2left_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_right2left_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_right2left_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_freq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
                    
                    n_infreq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2freq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2left_AL4_ind))
            
            elif sub  % 2 == 0:
                if bcount % 2 !=0:
                    
                    n_freq2infreq['1'].append(len(sw_right2left_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_right2left_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_right2left_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_left2right_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_left2right_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_left2right_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_left2right_AL4_ind))
                
                    n_freq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                   
                    n_infreq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
    
                    n_vertic2freq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2left_AL4_ind))
                    
                    
                    
    
                elif bcount % 2 ==0:
                    n_freq2infreq['1'].append(len(sw_left2right_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_left2right_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_left2right_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_left2right_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_right2left_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_right2left_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_right2left_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_freq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
                    
                    n_infreq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2freq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2left_AL4_ind))
                    
        elif sub > 8:
            if sub  % 2 != 0:
                if bcount % 2 !=0:
                    
                    n_freq2infreq['1'].append(len(sw_right2left_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_right2left_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_right2left_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_left2right_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_left2right_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_left2right_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_left2right_AL4_ind))
                
                    n_freq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                   
                    n_infreq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
    
                    n_vertic2freq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2left_AL4_ind))
                    
    
                elif bcount % 2 ==0:
                    n_freq2infreq['1'].append(len(sw_left2right_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_left2right_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_left2right_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_left2right_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_right2left_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_right2left_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_right2left_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_freq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
                    
                    n_infreq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2freq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2left_AL4_ind))
            
            elif sub  % 2 == 0:
                if bcount % 2 ==0:
                    
                    n_freq2infreq['1'].append(len(sw_right2left_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_right2left_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_right2left_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_left2right_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_left2right_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_left2right_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_left2right_AL4_ind))
                
                    n_freq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                   
                    n_infreq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
    
                    n_vertic2freq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2left_AL4_ind))
                    
                    
                    
    
                elif bcount % 2 !=0:
                    n_freq2infreq['1'].append(len(sw_left2right_AL1_ind))
                    n_freq2infreq['2'].append(len(sw_left2right_AL2_ind))
                    n_freq2infreq['3'].append(len(sw_left2right_AL3_ind))
                    n_freq2infreq['4'].append(len(sw_left2right_AL4_ind))
                    
                    n_infreq2freq['1'].append(len(sw_right2left_AL1_ind))
                    n_infreq2freq['2'].append(len(sw_right2left_AL2_ind))
                    n_infreq2freq['3'].append(len(sw_right2left_AL3_ind))
                    n_infreq2freq['4'].append(len(sw_right2left_AL4_ind))
                    
                    n_freq2vertic['1'].append(len(sw_left2vertic_AL1_ind))
                    n_freq2vertic['2'].append(len(sw_left2vertic_AL2_ind))
                    n_freq2vertic['3'].append(len(sw_left2vertic_AL3_ind))
                    n_freq2vertic['4'].append(len(sw_left2vertic_AL4_ind))
                    
                    n_infreq2vertic['1'].append(len(sw_right2vertic_AL1_ind))
                    n_infreq2vertic['2'].append(len(sw_right2vertic_AL2_ind))
                    n_infreq2vertic['3'].append(len(sw_right2vertic_AL3_ind))
                    n_infreq2vertic['4'].append(len(sw_right2vertic_AL4_ind))
                    
                    n_vertic2infreq['1'].append(len(sw_vertic2right_AL1_ind))
                    n_vertic2infreq['2'].append(len(sw_vertic2right_AL2_ind))
                    n_vertic2infreq['3'].append(len(sw_vertic2right_AL3_ind))
                    n_vertic2infreq['4'].append(len(sw_vertic2right_AL4_ind))
                    
                    n_vertic2freq['1'].append(len(sw_vertic2left_AL1_ind))
                    n_vertic2freq['2'].append(len(sw_vertic2left_AL2_ind))
                    n_vertic2freq['3'].append(len(sw_vertic2left_AL3_ind))
                    n_vertic2freq['4'].append(len(sw_vertic2left_AL4_ind))
                    

                
        
        '''indices for sensitivity analysis'''
        h_left_1_ind=all_indices((2,1,1),sa_zip)
        h_right_1_ind=all_indices((3,1,1),sa_zip)
        f1_ind=all_indices((1,0,1),sa_zip) 
        c1_ind=all_indices((1,1,1),sa_zip)
        m_left_1_ind=all_indices((2,0,1),sa_zip)
        m_right_1_ind=all_indices((3,0,1),sa_zip)

        
        h_left_2_ind=all_indices((2,1,2),sa_zip)
        h_right_2_ind=all_indices((3,1,2),sa_zip)
        f2_ind=all_indices((1,0,2),sa_zip) 
        c2_ind=all_indices((1,1,2),sa_zip)
        m_left_2_ind=all_indices((2,0,2),sa_zip)
        m_right_2_ind=all_indices((3,0,2),sa_zip)
        
        h_left_3_ind=all_indices((2,1,3),sa_zip)
        h_right_3_ind=all_indices((3,1,3),sa_zip)
        f3_ind=all_indices((1,0,3),sa_zip) 
        c3_ind=all_indices((1,1,3),sa_zip)
        m_left_3_ind=all_indices((2,0,3),sa_zip)
        m_right_3_ind=all_indices((3,0,3),sa_zip)
        
        h_left_4_ind=all_indices((2,1,4),sa_zip)
        h_right_4_ind=all_indices((3,1,4),sa_zip)
        f4_ind=all_indices((1,0,4),sa_zip) 
        c4_ind=all_indices((1,1,4),sa_zip)
        m_left_4_ind=all_indices((2,0,4),sa_zip)
        m_right_4_ind=all_indices((3,0,4),sa_zip)
        
        
        if block %2 ==0:
            SA_even['Hit1'].append(len(h_left_1_ind)+len(h_right_1_ind))
            SA_even['Hit2'].append(len(h_left_2_ind)+len(h_right_2_ind))  
            SA_even['Hit3'].append(len(h_left_3_ind)+len(h_right_3_ind))
            SA_even['Hit4'].append(len(h_left_4_ind)+len(h_right_4_ind))
            SA_even['Missed_Hit1'].append(len(m_left_1_ind)+len(m_right_1_ind))
            SA_even['Missed_Hit2'].append(len(m_left_2_ind)+len(m_right_2_ind))
            SA_even['Missed_Hit3'].append(len(m_left_3_ind)+len(m_right_3_ind))
            SA_even['Missed_Hit4'].append(len(m_left_4_ind)+len(m_right_4_ind))
            SA_even['False_Alarm1'].append(len(f1_ind))
            SA_even['False_Alarm2'].append(len(f2_ind))
            SA_even['False_Alarm3'].append(len(f3_ind))
            SA_even['False_Alarm4'].append(len(f4_ind))
            SA_even['Correct_Reject1'].append(len(c1_ind))
            SA_even['Correct_Reject2'].append(len(c2_ind))
            SA_even['Correct_Reject3'].append(len(c3_ind))
            SA_even['Correct_Reject4'].append(len(c4_ind))
        
        elif block % 2 != 0:
            SA_uneven['Hit1'].append(len(h_left_1_ind)+len(h_right_1_ind))
            SA_uneven['Hit2'].append(len(h_left_2_ind)+len(h_right_2_ind))
            SA_uneven['Hit3'].append(len(h_left_3_ind)+len(h_right_3_ind))
            SA_uneven['Hit4'].append(len(h_left_4_ind)+len(h_right_4_ind))
            SA_uneven['Missed_Hit1'].append(len(m_left_1_ind)+len(m_right_1_ind))
            SA_uneven['Missed_Hit2'].append(len(m_left_2_ind)+len(m_right_2_ind))
            SA_uneven['Missed_Hit3'].append(len(m_left_3_ind)+len(m_right_3_ind))
            SA_uneven['Missed_Hit4'].append(len(m_left_4_ind)+len(m_right_4_ind))
            SA_uneven['False_Alarm1'].append(len(f1_ind))
            SA_uneven['False_Alarm2'].append(len(f2_ind))
            SA_uneven['False_Alarm3'].append(len(f3_ind))
            SA_uneven['False_Alarm4'].append(len(f4_ind))
            SA_uneven['Correct_Reject1'].append(len(c1_ind))
            SA_uneven['Correct_Reject2'].append(len(c2_ind))
            SA_uneven['Correct_Reject3'].append(len(c3_ind))
            SA_uneven['Correct_Reject4'].append(len(c4_ind))
            

 
        '''no change trials'''
        condi_zip1=list(zip(ori,awareness_levels,switch))
        aware_zip1=list(zip(dummy,awareness_levels,switch))
        
        all_AL1_indices=all_indices((1,1),aware_zip)
        all_AL2_indices=all_indices((1,2),aware_zip)
        all_AL3_indices=all_indices((1,3),aware_zip)
        all_AL4_indices=all_indices((1,4),aware_zip)
        
        aware_lvls['all_AL1'].append(len(all_AL1_indices))
        aware_lvls['all_AL2'].append(len(all_AL2_indices))
        aware_lvls['all_AL3'].append(len(all_AL3_indices))
        aware_lvls['all_AL4'].append(len(all_AL4_indices))
        
        all_nan_recog_indices=all_indices((1,0),aware_zip)
        
        left_AL1_noswitch_indices=all_indices((2,1,0),condi_zip1)
        right_AL1_noswitch_indices=all_indices((3,1,0),condi_zip1)
        vertic_AL1_noswitch_indices=all_indices((1,1,0),condi_zip1) 
       
        left_AL2_noswitch_indices=all_indices((2,2,0),condi_zip1)
        right_AL2_noswitch_indices=all_indices((3,2,0),condi_zip1)
        vertic_AL2_noswitch_indices=all_indices((1,2,0),condi_zip1) 
        
        left_AL3_noswitch_indices=all_indices((2,3,0),condi_zip1)
        right_AL3_noswitch_indices=all_indices((3,3,0),condi_zip1)
        vertic_AL3_noswitch_indices=all_indices((1,3,0),condi_zip1) 
        
        left_AL4_noswitch_indices=all_indices((2,4,0),condi_zip1)
        right_AL4_noswitch_indices=all_indices((3,4,0),condi_zip1)
        vertic_AL4_noswitch_indices=all_indices((1,4,0),condi_zip1) 
        
        left_AL_nan=all_indices((2,0),condi_zip1)
        right_AL_nan=all_indices((3,0),condi_zip1)
        vertic_AL_nan=all_indices((1,0),condi_zip1) 
        
        left_AL_nan=all_indices((2,1),condi_zip1)
        right_AL_nan=all_indices((3,1),condi_zip1)
        vertic_AL_nan=all_indices((1,1),condi_zip1) 
        
   
# =====================================================================================
#               Get the values for each condition 
#               (assign values to the keys of all dics by using my TrialFinder Function)
# ====================================================================================
        '''switch depending on prior ori'''
        if sub < 8:
            if sub % 2 ==0: # even blocks are right-weighted!
                
                '''collect n of trials'''

                '''switch'''
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,freq2infreq_rt['1'],sub_rt,block)#takes the first ind on even blocks and the second ind list on uneven blocks
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,freq2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,freq2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,freq2infreq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,infreq2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,infreq2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,infreq2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,infreq2freq_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,freq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,freq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,freq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,freq2vertic_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,infreq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,infreq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,infreq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,infreq2vertic_rt['4'],sub_rt,block)
                
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2freq_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2infreq_rt['4'],sub_rt,block)
                
                
                '''no switch'''
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_freq_rt['1'],sub_rt,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_freq_rt['2'],sub_rt,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_freq_rt['3'],sub_rt,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_freq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_infreq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_rt['4'],sub_rt,block)
            
            
            
            
            elif sub % 2 !=0: #even blocks are left-weighted
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,freq2infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,freq2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,freq2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,freq2infreq_rt['4'],sub_rt,block)
        
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,infreq2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,infreq2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,infreq2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,infreq2freq_rt['4'],sub_rt,block)
                    
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,freq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,freq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,freq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,freq2vertic_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,infreq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,infreq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,infreq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,infreq2vertic_rt['4'],sub_rt,block)
                
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2freq_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2infreq_rt['4'],sub_rt,block)
                
                '''no switch'''
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_freq_rt['1'],sub_rt,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_freq_rt['2'],sub_rt,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_freq_rt['3'],sub_rt,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_freq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_infreq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_rt['4'],sub_rt,block)
            
        
        elif sub >= 8:
            if sub % 2 !=0: # even blocks are right-weighted!
                
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,freq2infreq_rt['1'],sub_rt,block)#even block right if sub uneven; changed with 9th subject! 
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,freq2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,freq2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,freq2infreq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,infreq2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,infreq2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,infreq2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,infreq2freq_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,freq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,freq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,freq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,freq2vertic_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,infreq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,infreq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,infreq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,infreq2vertic_rt['4'],sub_rt,block)
                
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2freq_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2infreq_rt['4'],sub_rt,block)

                '''no switch'''
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_freq_rt['1'],sub_rt,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_freq_rt['2'],sub_rt,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_freq_rt['3'],sub_rt,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_freq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_infreq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_rt['4'],sub_rt,block)
            
                
            
            
            elif sub % 2 ==0: 
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,freq2infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,freq2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,freq2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,freq2infreq_rt['4'],sub_rt,block)
        
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,infreq2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,infreq2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,infreq2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,infreq2freq_rt['4'],sub_rt,block)
                    
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,freq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,freq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,freq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,freq2vertic_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,infreq2vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,infreq2vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,infreq2vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,infreq2vertic_rt['4'],sub_rt,block)
                
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2freq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2freq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2freq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2freq_rt['4'],sub_rt,block)
            
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2infreq_rt['4'],sub_rt,block)
           
                                              
                '''no switch'''
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_freq_rt['1'],sub_rt,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_freq_rt['2'],sub_rt,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_freq_rt['3'],sub_rt,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_freq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_infreq_rt['1'],sub_rt,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_infreq_rt['2'],sub_rt,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_infreq_rt['3'],sub_rt,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_infreq_rt['4'],sub_rt,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_rt['1'],sub_rt,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_rt['2'],sub_rt,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_rt['3'],sub_rt,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_rt['4'],sub_rt,block)


    '''flatten all rt lists and remove zeros (trimmed rts) - these lists count the number of trials! '''
    
    '''with 0.2 I get rid of any reaction time slower than 200ms - some subjects might initiate their response while the mask is still 
    on screen but it is collected in the very beginning of the actual response window leading to extremly short (invalid) reaction times'''
    

    for key in freq2infreq_rt:
        freq2infreq_rt[key]=[item for sublist in freq2infreq_rt[key] for item in sublist] #flattening, removes empty lists
        freq2infreq_rt[key]=[x for x in freq2infreq_rt[key] if x >0.2] #remove nans/ RTs < 200ms
        
   
    for key in infreq2freq_rt:
        infreq2freq_rt[key]=[item for sublist in infreq2freq_rt[key] for item in sublist] 
        infreq2freq_rt[key]=[x for x in infreq2freq_rt[key] if x >0.2]

        
    for key in freq2vertic_rt:
        freq2vertic_rt[key]=[item for sublist in freq2vertic_rt[key] for item in sublist] 
        freq2vertic_rt[key]=[x for x in freq2vertic_rt[key] if x >0.2] 

    
    for key in infreq2vertic_rt:
        infreq2vertic_rt[key]=[item for sublist in infreq2vertic_rt[key] for item in sublist] 
        infreq2vertic_rt[key]=[x for x in infreq2vertic_rt[key] if x >0.2] 

        
    for key in vertic2infreq_rt:
        vertic2infreq_rt[key]=[item for sublist in vertic2infreq_rt[key] for item in sublist] 
        vertic2infreq_rt[key]=[x for x in vertic2infreq_rt[key] if x >0.2] 

        
        
    for key in vertic2freq_rt:
        vertic2freq_rt[key]=[item for sublist in vertic2freq_rt[key] for item in sublist] 
        vertic2freq_rt[key]=[x for x in vertic2freq_rt[key] if x >0.2] 

       
 
    '''NOTE: now the script is changed so that I allow for responses with Mask onset (actual reaction times) '''
    '''use dictionaries of no switch trials'''
    
    for key in ns_freq_rt:
        ns_freq_rt[key]=[item for sublist in ns_freq_rt[key] for item in sublist] #flattening
        ns_freq_rt[key]=[x for x in ns_freq_rt[key] if x >0.2] #remove nans


    for key in ns_infreq_rt:
        ns_infreq_rt[key]=[item for sublist in ns_infreq_rt[key] for item in sublist] #flattening
        ns_infreq_rt[key]=[x for x in ns_infreq_rt[key] if x >0.2] #remove nans

        
    for key in ns_vertic_rt:
        ns_vertic_rt[key]=[item for sublist in ns_vertic_rt[key] for item in sublist] #flattening
        ns_vertic_rt[key]=[x for x in ns_vertic_rt[key] if x >0.2] #remove nans

       
   
    print ('All RT lists were flattened')
    
    '''sum up number of trials '''
    for key in n_freq2infreq:
        n_freq2infreq[key]=sum(n_freq2infreq[key])
    
        
        
        
        
# =============================================================================
#     Calculations 
# =============================================================================
    


    '''calculate mean RTs, here work around empty lists'''

    #switch
    for key in freq2infreq_rt:
        if len(freq2infreq_rt[key])> 1:freq2infreq_rt[key]=trimm(freq2infreq_rt[key])
        else:freq2infreq_rt[key]=np.nanmean(freq2infreq_rt[key])
    mean_freq2infreq_rt=freq2infreq_rt
    
    for key in infreq2freq_rt:
        if len(infreq2freq_rt[key])> 1:
            infreq2freq_rt[key]=trimm(infreq2freq_rt[key])
        else:infreq2freq_rt[key]=np.nanmean(infreq2freq_rt[key])
    mean_infreq2freq_rt=infreq2freq_rt
     

    for key in infreq2vertic_rt:
        if len(infreq2vertic_rt[key])>1:
            infreq2vertic_rt[key]=trimm(infreq2vertic_rt[key])#
        else:infreq2vertic_rt[key]=np.nanmean(infreq2vertic_rt[key])
    mean_infreq2vertic_rt=infreq2vertic_rt
    
    
    for key in freq2vertic_rt:
        if len(freq2vertic_rt[key])>1:
            freq2vertic_rt[key]=trimm(freq2vertic_rt[key])#
        else:freq2vertic_rt[key]=np.nanmean(freq2vertic_rt[key])
    mean_freq2vertic_rt=freq2vertic_rt
    
    for key in vertic2infreq_rt:
        if len(vertic2infreq_rt[key])>1:
            vertic2infreq_rt[key]=trimm(vertic2infreq_rt[key])#
        else:vertic2infreq_rt[key]=np.nanmean(vertic2infreq_rt[key])
    mean_vertic2infreq_rt=vertic2infreq_rt
    
    for key in vertic2freq_rt:
        if len(vertic2freq_rt[key])>1:
            vertic2freq_rt[key]=trimm(vertic2freq_rt[key])#
        else:vertic2freq_rt[key]=np.nanmean(vertic2freq_rt[key])
    mean_vertic2freq_rt=vertic2freq_rt
   
    
    
    #no switch
    for key in ns_vertic_rt:
        if len(ns_vertic_rt[key])> 1:
            ns_vertic_rt[key]=trimm(ns_vertic_rt[key])#
        else:ns_vertic_rt[key]=np.mean(ns_vertic_rt[key])
    mean_ns_vertic_rt=ns_vertic_rt

    for key in ns_freq_rt:
        if len(ns_freq_rt[key])> 1:
            ns_freq_rt[key]=trimm(ns_freq_rt[key])#
        else:ns_freq_rt[key]=np.nanmean(ns_freq_rt[key])
    mean_ns_freq_rt=ns_freq_rt
    
    
    for key in ns_infreq_rt:
        if len(ns_infreq_rt[key])> 1:
            ns_infreq_rt[key]=trimm(ns_infreq_rt[key])#
        else:ns_infreq_rt[key]=np.nanmean(ns_infreq_rt[key])
    mean_ns_infreq_rt=ns_infreq_rt
    

    
    '''WRITE CSV FILE CONTAINING RTs (order is important here!)'''
    
    params=['sub','rt_freq2infreq_1','rt_freq2infreq_2','rt_freq2infreq_3','rt_freq2infreq_4',
              'rt_freq2vertic_1','rt_freq2vertic_2','rt_freq2vertic_3','rt_freq2vertic_4',
              
              'rt_infreq2freq_1','rt_infreq2freq_2','rt_infreq2freq_3','rt_infreq2freq_4',
              'rt_infreq2vertic_1','rt_infreq2vertic_2','rt_infreq2vertic_3','rt_infreq2vertic_4',
              
              'rt_vertic2freq_1','rt_ivertic2freq_2','rt_vertic2freq_3','rt_vertic2freq_4',
              'rt_vertic2infreq_1','rt_vertic2infreq_2','rt_vertic2infreq_3','rt_vertic2infreq_4',
              
              'rt_freq_ns_1','rt_freq_ns_2','rt_freq_ns_3','rt_freq_ns_4',
               'rt_infreq_ns_1','rt_infreq_ns_2','rt_infreq_ns_3','rt_infreq_ns_4',
               'rt_vertic_ns_1','rt_vertic_ns_2','rt_vertic_ns_3','rt_vertic_ns_4']

             
              
    rts=[sub,mean_freq2infreq_rt['1'],mean_freq2infreq_rt['2'],mean_freq2infreq_rt['3'],mean_freq2infreq_rt['4'],
                mean_freq2vertic_rt['1'],mean_freq2vertic_rt['2'],mean_freq2vertic_rt['3'],mean_freq2vertic_rt['4'],
                mean_infreq2freq_rt['1'],mean_infreq2freq_rt['2'],mean_infreq2freq_rt['3'],mean_infreq2freq_rt['4'],
                mean_infreq2vertic_rt['1'],mean_infreq2vertic_rt['2'],mean_infreq2vertic_rt['3'],mean_infreq2vertic_rt['4'],
                mean_vertic2freq_rt['1'],mean_vertic2freq_rt['2'],mean_vertic2freq_rt['3'],mean_vertic2freq_rt['4'],
                mean_vertic2infreq_rt['1'],mean_vertic2infreq_rt['2'],mean_vertic2infreq_rt['3'],mean_vertic2infreq_rt['4'],
                mean_ns_freq_rt['1'],mean_ns_freq_rt['2'],mean_ns_freq_rt['3'],mean_ns_freq_rt['4'],
                mean_ns_infreq_rt['1'],mean_ns_infreq_rt['2'],mean_ns_infreq_rt['3'],mean_ns_infreq_rt['4'],
                mean_ns_vertic_rt['1'],mean_ns_vertic_rt['2'],mean_ns_vertic_rt['3'],mean_ns_vertic_rt['4']
                ] 
#                   
###    
    
    
    rt_spec_sw_al1=np.nanmean([mean_freq2infreq_rt['1'],mean_freq2vertic_rt['1']])
    rt_av_sw_al1=np.nanmean([mean_freq2infreq_rt['1'],mean_freq2vertic_rt['1'],mean_infreq2freq_rt['1'],mean_infreq2vertic_rt['1'],
                             mean_vertic2freq_rt['1'],mean_vertic2infreq_rt['1']])
    
    
   
    RT.append(rt_spec_sw_al1)
    avRT.append(rt_av_sw_al1)
    AL.append('AL1')
    SW.append('switch')
    SUB.append(sub)
    
    rt_spec_sw_al2=np.nanmean([mean_freq2infreq_rt['2'],mean_freq2vertic_rt['2']])
    rt_av_sw_al2=np.nanmean([mean_freq2infreq_rt['2'],mean_freq2vertic_rt['2'],mean_infreq2freq_rt['2'],mean_infreq2vertic_rt['2'],
                             mean_vertic2freq_rt['2'],mean_vertic2infreq_rt['2']])

    RT.append(rt_spec_sw_al2)
    avRT.append(rt_av_sw_al2)
    AL.append('AL2')
    SW.append('switch')
    SUB.append(sub)
    
    rt_spec_sw_al3=np.nanmean([mean_freq2infreq_rt['3'],mean_freq2vertic_rt['3']])
    rt_av_sw_al3=np.nanmean([mean_freq2infreq_rt['3'],mean_freq2vertic_rt['3'],mean_infreq2freq_rt['3'],mean_infreq2vertic_rt['3'],
                             mean_vertic2freq_rt['3'],mean_vertic2infreq_rt['3']])

    RT.append(rt_spec_sw_al3)
    avRT.append(rt_av_sw_al3)
    AL.append('AL3')
    SW.append('switch')
    SUB.append(sub)
    
    #repeat
    ns_al1=np.nanmean([mean_ns_vertic_rt['1'],mean_ns_infreq_rt['1'],mean_ns_freq_rt['1']])
    RT.append(ns_al1)
    avRT.append(ns_al1)
    AL.append('AL1')
    SW.append('repeat')
    SUB.append(sub)
    ns_al2=np.nanmean([mean_ns_vertic_rt['2'],mean_ns_infreq_rt['2'] ,mean_ns_freq_rt['2']])
    RT.append(ns_al2)
    avRT.append(ns_al2)
    AL.append('AL2')
    SW.append('repeat')
    SUB.append(sub)
    ns_al3=np.nanmean([mean_ns_vertic_rt['3'], mean_ns_infreq_rt['3'] , mean_ns_freq_rt['3']])
    RT.append(ns_al3)
    avRT.append(ns_al3)
    AL.append('AL3')
    SW.append('repeat')
    SUB.append(sub)
    
    
    
    #sw costs
    
    s_cost1=rt_spec_sw_al1 - ns_al1
    s_cost2=rt_spec_sw_al2 - ns_al2
    s_cost3=rt_spec_sw_al3 - ns_al3
    
    a_cost1=rt_av_sw_al1 - ns_al1
    a_cost2=rt_av_sw_al2 - ns_al2
    a_cost3=rt_av_sw_al3 - ns_al3
    
    lb=['sub','spec_costs1','spec_costs2','spec_costs3','av_costs1','av_costs2','av_costs3']
    costs=[sub,s_cost3,s_cost3,s_cost3,a_cost1,a_cost2,a_cost3]
    

    df=list(zip(SUB,RT,AL,SW,avRT))
    
    df=pd.DataFrame(df)
    now = datetime.now()
    date = now.strftime("%m_%d_%Y")
    df.to_csv(datapath1 + 'pretarget_mean_RTs_long_allsubs_%s.csv'%(date), sep=',', index= False, header=True)


 