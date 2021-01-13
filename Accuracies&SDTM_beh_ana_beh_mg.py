#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:32:22 2019

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

datapath1=wdir+ u'/results/Accuracies_and_SDTM/'

'''if sub is even: block is even = right-weighted, block = uneven left_weighted
if sub is uneven: block is even = left_weightes, block = uneven right_weighted'''  

nTrials=36

count_nan=0


'''
sub1 gerade blöcke -> 3 also rechts, ungerade 2 also links '''

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
    # Floors an ceilings are replaced by half hits and half FA's
    try:
        halfHit = 0.5/(hits+misses)
    except ZeroDivisionError:
         halfHit = 0.5
        
    try:
        halfFa = 0.5/(fas+crs)
    except ZeroDivisionError:
        halfFa = 0.5
 
    # Calculate hitrate and avoid d'infinity
    try:
        hitRate = hits/(hits+misses)
    except ZeroDivisionError:
        hitRate = 0
        
    if hitRate == 1: hitRate = 1-halfHit
    if hitRate == 0: hitRate = halfHit
 
    # Calculate false alarm rate and avoid d' infinity
    try: 
        faRate = fas/(fas+crs)
    except ZeroDivisionError:
        faRate= 0
            
    if faRate == 1: faRate = 1-halfFa
    if faRate == 0: faRate = halfFa
 
    # Return d', beta, c , P(c), Ad' and A  # these are the formulas for yes/no task sensitivity (see Macmillan)
    # c = criterion location used as a measure of bias
    # A = our sensitivity measure
    #p(c) = 
    out = {}
    out['d'] = Z(hitRate) - Z(faRate) #that´s for a yes/no task!!!
    out['beta'] = exp((Z(faRate)**2 - Z(hitRate)**2)/2)
    out['c'] = -(Z(hitRate) + Z(faRate))/2
    out['Ad'] = norm.cdf(out['d']/sqrt(2))

    try:
        out['Pc']=0.444*hitRate+0.222*hitRate+0.333*(1-faRate)
    except ZeroDivisionError:
        out['Pc']='nan'
    if hitRate >= faRate:
        out['A']=0.5+(((hitRate-faRate)*(1+hitRate-faRate))/(4*hitRate+(1-faRate)))
    elif hitRate < faRate:
        out['A']=0.5-(((faRate-hitRate)*(1+faRate-hitRate))/(4*faRate+(1-hitRate)))
        
    return out 
    #print (out)
    

ACC=[]
avACC=[]
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
    freq2infreq_acc={}
    infreq2freq_acc={}
    freq2vertic_acc={}
    infreq2vertic_acc={}
    vertic2freq_acc={}
    vertic2infreq_acc={}
    infreq2vertic_acc={}
    freq2vertic_acc={}
    ns_freq_acc={}
    ns_infreq_acc={}
    ns_vertic_acc={}
       
    ''' fill the dictionary, every entry becomes an empty list for data collection. 
    These can be flattened and used for mean calculation. '''
    
    lvls=[1,2,3,4]
    for l in lvls:
        freq2infreq_acc['{0}'.format(l)]=[]
        infreq2freq_acc['{0}'.format(l)]=[]
        freq2vertic_acc['{0}'.format(l)]=[]
        infreq2vertic_acc['{0}'.format(l)]=[]
        vertic2freq_acc['{0}'.format(l)]=[]
        vertic2infreq_acc['{0}'.format(l)]=[]
        ns_freq_acc['{0}'.format(l)]=[]
        ns_infreq_acc['{0}'.format(l)]=[]
        ns_vertic_acc['{0}'.format(l)]=[]
        infreq2vertic_acc['{0}'.format(l)]=[]
        freq2vertic_acc['{0}'.format(l)]=[]
        

    
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
        left=0
        
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
                left=left+1
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
                
        #print (sub_acc)  
        print ('sub',sub) 
        print ('block',block)
        print ('left',left)
        
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
        
        
        # sensitivity analysis
        h_left_1_ind=all_indices((2,1,1),sa_zip)#hit 1=correct left-tilt 
        h_right_1_ind=all_indices((3,1,1),sa_zip)#hit 2= correct right-tilt
        f1_ind=all_indices((1,0,1),sa_zip)#FA=incorrect vertical 
        c1_ind=all_indices((1,1,1),sa_zip)#CR=correct vertical
        m_left_1_ind=all_indices((2,0,1),sa_zip)#MH1= incorrect tilt left
        m_right_1_ind=all_indices((3,0,1),sa_zip)#MH2=incorrect tilt right

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
    
   
# =====================================================================================
#               Get the values for each condition 
#               (assign values to the keys of all dics by using my TrialFinder Function)
# ====================================================================================
        '''switch depending on prior ori'''
        if sub < 8:
            if sub % 2 ==0: # even blocks are right-weighted!
                
                '''collect accuracies for each condition (type of switch and repeat) blockwise'''

                '''switch'''
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,freq2infreq_acc['1'],sub_acc,block)#takes the first ind on even blocks and the second ind list on uneven blocks
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,freq2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,freq2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,freq2infreq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,infreq2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,infreq2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,infreq2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,infreq2freq_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,freq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,freq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,freq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,freq2vertic_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,infreq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,infreq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,infreq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,infreq2vertic_acc['4'],sub_acc,block)
                
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2freq_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2infreq_acc['4'],sub_acc,block)
                
                
                '''no switch'''
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_freq_acc['1'],sub_acc,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_freq_acc['2'],sub_acc,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_freq_acc['3'],sub_acc,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_freq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_infreq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_acc['4'],sub_acc,block)
            
            
            
            
            elif sub % 2 !=0: #even blocks are left-weighted
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,freq2infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,freq2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,freq2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,freq2infreq_acc['4'],sub_acc,block)
        
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,infreq2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,infreq2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,infreq2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,infreq2freq_acc['4'],sub_acc,block)
                    
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,freq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,freq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,freq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,freq2vertic_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,infreq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,infreq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,infreq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,infreq2vertic_acc['4'],sub_acc,block)
                
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2freq_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2infreq_acc['4'],sub_acc,block)
                
                '''no switch'''
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_freq_acc['1'],sub_acc,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_freq_acc['2'],sub_acc,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_freq_acc['3'],sub_acc,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_freq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_infreq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_acc['4'],sub_acc,block)
            
        
        elif sub >= 8:
            if sub % 2 !=0: # even blocks are right-weighted!
                
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,freq2infreq_acc['1'],sub_acc,block)#even block right if sub uneven; changed with 9th subject! 
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,freq2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,freq2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,freq2infreq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,infreq2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,infreq2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,infreq2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,infreq2freq_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,freq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,freq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,freq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,freq2vertic_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,infreq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,infreq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,infreq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,infreq2vertic_acc['4'],sub_acc,block)
                
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2freq_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2infreq_acc['4'],sub_acc,block)

                '''no switch'''
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_freq_acc['1'],sub_acc,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_freq_acc['2'],sub_acc,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_freq_acc['3'],sub_acc,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_freq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_infreq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_acc['4'],sub_acc,block)
            
                
            
            
            elif sub % 2 ==0: 
                
                TrialFinder_switch(sw_left2right_AL1_ind,sw_right2left_AL1_ind,freq2infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL2_ind,sw_right2left_AL2_ind,freq2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL3_ind,sw_right2left_AL3_ind,freq2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2right_AL4_ind,sw_right2left_AL4_ind,freq2infreq_acc['4'],sub_acc,block)
        
                TrialFinder_switch(sw_right2left_AL1_ind,sw_left2right_AL1_ind,infreq2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL2_ind,sw_left2right_AL2_ind,infreq2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL3_ind,sw_left2right_AL3_ind,infreq2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2left_AL4_ind,sw_left2right_AL4_ind,infreq2freq_acc['4'],sub_acc,block)
                    
                TrialFinder_switch(sw_left2vertic_AL1_ind,sw_right2vertic_AL1_ind,freq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL2_ind,sw_right2vertic_AL2_ind,freq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL3_ind,sw_right2vertic_AL3_ind,freq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_left2vertic_AL4_ind,sw_right2vertic_AL4_ind,freq2vertic_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_right2vertic_AL1_ind,sw_left2vertic_AL1_ind,infreq2vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL2_ind,sw_left2vertic_AL2_ind,infreq2vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL3_ind,sw_left2vertic_AL3_ind,infreq2vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_right2vertic_AL4_ind,sw_left2vertic_AL4_ind,infreq2vertic_acc['4'],sub_acc,block)
                
                TrialFinder_switch(sw_vertic2left_AL1_ind,sw_vertic2right_AL1_ind,vertic2freq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL2_ind,sw_vertic2right_AL2_ind,vertic2freq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL3_ind,sw_vertic2right_AL3_ind,vertic2freq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2left_AL4_ind,sw_vertic2right_AL4_ind,vertic2freq_acc['4'],sub_acc,block)
            
                TrialFinder_switch(sw_vertic2right_AL1_ind,sw_vertic2left_AL1_ind,vertic2infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL2_ind,sw_vertic2left_AL2_ind,vertic2infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL3_ind,sw_vertic2left_AL3_ind,vertic2infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(sw_vertic2right_AL4_ind,sw_vertic2left_AL4_ind,vertic2infreq_acc['4'],sub_acc,block)
           
                                              
                '''no switch'''
                TrialFinder_switch(left_AL1_noswitch_indices,right_AL1_noswitch_indices,ns_freq_acc['1'],sub_acc,block)
                TrialFinder_switch(left_AL2_noswitch_indices,right_AL2_noswitch_indices,ns_freq_acc['2'],sub_acc,block)
                TrialFinder_switch(left_AL3_noswitch_indices,right_AL3_noswitch_indices,ns_freq_acc['3'],sub_acc,block)
                TrialFinder_switch(left_AL4_noswitch_indices,right_AL4_noswitch_indices,ns_freq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(right_AL1_noswitch_indices,left_AL1_noswitch_indices,ns_infreq_acc['1'],sub_acc,block)
                TrialFinder_switch(right_AL2_noswitch_indices,left_AL2_noswitch_indices,ns_infreq_acc['2'],sub_acc,block)
                TrialFinder_switch(right_AL3_noswitch_indices,left_AL3_noswitch_indices,ns_infreq_acc['3'],sub_acc,block)
                TrialFinder_switch(right_AL4_noswitch_indices,left_AL4_noswitch_indices,ns_infreq_acc['4'],sub_acc,block)
                
                TrialFinder_switch(vertic_AL1_noswitch_indices,vertic_AL1_noswitch_indices,ns_vertic_acc['1'],sub_acc,block)
                TrialFinder_switch(vertic_AL2_noswitch_indices,vertic_AL2_noswitch_indices,ns_vertic_acc['2'],sub_acc,block)
                TrialFinder_switch(vertic_AL3_noswitch_indices,vertic_AL3_noswitch_indices,ns_vertic_acc['3'],sub_acc,block)
                TrialFinder_switch(vertic_AL4_noswitch_indices,vertic_AL4_noswitch_indices,ns_vertic_acc['4'],sub_acc,block)


    '''flatten all acc lists to get one list containing all accs of all blocks '''
    
    
    for key in freq2infreq_acc:
        freq2infreq_acc[key]=[item for sublist in freq2infreq_acc[key] for item in sublist] #flattening, removes empty lists
    

   
    for key in infreq2freq_acc:
        infreq2freq_acc[key]=[item for sublist in infreq2freq_acc[key] for item in sublist] 
            
    for key in freq2vertic_acc:
        freq2vertic_acc[key]=[item for sublist in freq2vertic_acc[key] for item in sublist] 
         
    
    for key in infreq2vertic_acc:
        infreq2vertic_acc[key]=[item for sublist in infreq2vertic_acc[key] for item in sublist] 
              
    for key in vertic2infreq_acc:
        vertic2infreq_acc[key]=[item for sublist in vertic2infreq_acc[key] for item in sublist] 
             
        
    for key in vertic2freq_acc:
        vertic2freq_acc[key]=[item for sublist in vertic2freq_acc[key] for item in sublist] 
          
 
    '''use dictionaries of no switch trials'''
    
    for key in ns_freq_acc:
        ns_freq_acc[key]=[item for sublist in ns_freq_acc[key] for item in sublist] #flattening
               #uneven blocks same number of key
    for key in ns_infreq_acc:
        ns_infreq_acc[key]=[item for sublist in ns_infreq_acc[key] for item in sublist] #flattening
             
    for key in ns_vertic_acc:
        ns_vertic_acc[key]=[item for sublist in ns_vertic_acc[key] for item in sublist] #flattening
         
   
    print ('All RT lists were flattened')
    
         
    
        
        
        
# =============================================================================
#     Calculations 
# =============================================================================
    
      
    
    '''calculate mean ACCs'''

    #switch
    for key in freq2infreq_acc:
        if len(freq2infreq_acc[key])>0:
            freq2infreq_acc[key]=np.nanmean(freq2infreq_acc[key])
        else:freq2infreq_acc[key]='nan'
    mean_freq2infreq_acc=freq2infreq_acc
    
    for key in infreq2freq_acc:
        if len(infreq2freq_acc[key])> 0:
            infreq2freq_acc[key]=np.nanmean(infreq2freq_acc[key])
        else:infreq2freq_acc[key]='nan'
    mean_infreq2freq_acc=infreq2freq_acc
     

    for key in infreq2vertic_acc:
        if len(infreq2vertic_acc[key])>0:
            infreq2vertic_acc[key]=np.nanmean(infreq2vertic_acc[key])
        else:infreq2vertic_acc[key]='nan'
    mean_infreq2vertic_acc=infreq2vertic_acc
    
    
    for key in freq2vertic_acc:
        if len(freq2vertic_acc[key])>0:
            freq2vertic_acc[key]=np.nanmean(freq2vertic_acc[key])
        else:freq2vertic_acc[key]='nan'
    mean_freq2vertic_acc=freq2vertic_acc
    
    for key in vertic2infreq_acc:
        if len(vertic2infreq_acc[key])>0:
            vertic2infreq_acc[key]=np.nanmean(vertic2infreq_acc[key])
        else:vertic2infreq_acc[key]='nan'
    mean_vertic2infreq_acc=vertic2infreq_acc
    
    for key in vertic2freq_acc:
        if len(vertic2freq_acc[key])>0:
            vertic2freq_acc[key]=np.nanmean(vertic2freq_acc[key])
        else:vertic2freq_acc[key]='nan'
    mean_vertic2freq_acc=vertic2freq_acc
   
    
    
    #no switch
    for key in ns_vertic_acc:
        if len(ns_vertic_acc[key])> 0:
            ns_vertic_acc[key]=np.mean(ns_vertic_acc[key])#
        else:ns_vertic_acc[key]='nan'
    mean_ns_vertic_acc=ns_vertic_acc

    for key in ns_freq_acc:
        if len(ns_freq_acc[key])> 0:
            ns_freq_acc[key]=np.nanmean(ns_freq_acc[key])
        else:ns_freq_acc[key]='nan'
    mean_ns_freq_acc=ns_freq_acc
    
    
    for key in ns_infreq_acc:
        if len(ns_infreq_acc[key])> 0:
            ns_infreq_acc[key]=np.nanmean(ns_infreq_acc[key])
        else:ns_infreq_acc[key]='nan'
    mean_ns_infreq_acc=ns_infreq_acc
    

    
    '''WRITE CSV FILE CONTAINING RTs (order is important here!)'''
    
    params=['sub','acc_freq2infreq_1','acc_freq2infreq_2','acc_freq2infreq_3','acc_freq2infreq_4',
              'acc_freq2veaccic_1','acc_freq2veaccic_2','acc_freq2veaccic_3','acc_freq2veaccic_4',
              
              'acc_infreq2freq_1','acc_infreq2freq_2','acc_infreq2freq_3','acc_infreq2freq_4',
              'acc_infreq2veaccic_1','acc_infreq2veaccic_2','acc_infreq2veaccic_3','acc_infreq2veaccic_4',
              
              'acc_veaccic2freq_1','acc_iveaccic2freq_2','acc_veaccic2freq_3','acc_veaccic2freq_4',
              'acc_veaccic2infreq_1','acc_veaccic2infreq_2','acc_veaccic2infreq_3','acc_veaccic2infreq_4',
              
              'acc_freq_ns_1','acc_freq_ns_2','acc_freq_ns_3','acc_freq_ns_4',
               'acc_infreq_ns_1','acc_infreq_ns_2','acc_infreq_ns_3','acc_infreq_ns_4',
               'acc_veaccic_ns_1','acc_veaccic_ns_2','acc_veaccic_ns_3','acc_veaccic_ns_4']

             
              
    rts=[sub,mean_freq2infreq_acc['1'],mean_freq2infreq_acc['2'],mean_freq2infreq_acc['3'],mean_freq2infreq_acc['4'],
                mean_freq2vertic_acc['1'],mean_freq2vertic_acc['2'],mean_freq2vertic_acc['3'],mean_freq2vertic_acc['4'],
                mean_infreq2freq_acc['1'],mean_infreq2freq_acc['2'],mean_infreq2freq_acc['3'],mean_infreq2freq_acc['4'],
                mean_infreq2vertic_acc['1'],mean_infreq2vertic_acc['2'],mean_infreq2vertic_acc['3'],mean_infreq2vertic_acc['4'],
                mean_vertic2freq_acc['1'],mean_vertic2freq_acc['2'],mean_vertic2freq_acc['3'],mean_vertic2freq_acc['4'],
                mean_vertic2infreq_acc['1'],mean_vertic2infreq_acc['2'],mean_vertic2infreq_acc['3'],mean_vertic2infreq_acc['4'],
                mean_ns_freq_acc['1'],mean_ns_freq_acc['2'],mean_ns_freq_acc['3'],mean_ns_freq_acc['4'],
                mean_ns_infreq_acc['1'],mean_ns_infreq_acc['2'],mean_ns_infreq_acc['3'],mean_ns_infreq_acc['4'],
                mean_ns_vertic_acc['1'],mean_ns_vertic_acc['2'],mean_ns_vertic_acc['3'],mean_ns_vertic_acc['4']
                ] 
#                   
###    
    
    
    acc_spec_sw_al1=np.array([mean_freq2infreq_acc['1'],mean_freq2vertic_acc['1']],np.float)
    acc_spec_sw_al1=np.nanmean(acc_spec_sw_al1)
    acc_av_sw_al1=np.array([mean_freq2infreq_acc['1'],mean_freq2vertic_acc['1'],mean_infreq2freq_acc['1'],mean_infreq2vertic_acc['1'],
                             mean_vertic2freq_acc['1'],mean_vertic2infreq_acc['1']],np.float)
    acc_av_sw_al1=np.nanmean(acc_av_sw_al1)
    
    
   
    
    ACC.append(acc_spec_sw_al1)
    avACC.append(acc_av_sw_al1)
    AL.append('AL1')
    SW.append('switch')
    SUB.append(sub)
    
    acc_spec_sw_al2=np.array([mean_freq2infreq_acc['2'],mean_freq2vertic_acc['2']],np.float)
    acc_spec_sw_al2=np.nanmean(acc_spec_sw_al2)
    acc_av_sw_al2=np.array([mean_freq2infreq_acc['2'],mean_freq2vertic_acc['2'],mean_infreq2freq_acc['2'],mean_infreq2vertic_acc['2'],
                             mean_vertic2freq_acc['2'],mean_vertic2infreq_acc['2']],np.float)
    acc_av_sw_al2=np.nanmean(acc_av_sw_al2)
    

    ACC.append(acc_spec_sw_al2)
    avACC.append(acc_av_sw_al2)
    AL.append('AL2')
    SW.append('switch')
    SUB.append(sub)
    
    acc_spec_sw_al3=np.array([mean_freq2infreq_acc['3'],mean_freq2vertic_acc['3']],np.float)
    acc_spec_sw_al3=np.nanmean(acc_spec_sw_al3)
    acc_av_sw_al3=np.array([mean_freq2infreq_acc['3'],mean_freq2vertic_acc['3'],mean_infreq2freq_acc['3'],mean_infreq2vertic_acc['3'],
                             mean_vertic2freq_acc['3'],mean_vertic2infreq_acc['3']],np.float)
    acc_av_sw_al3=np.nanmean(acc_av_sw_al3)
    

    ACC.append(acc_spec_sw_al3)
    avACC.append(acc_av_sw_al3)
    AL.append('AL3')
    SW.append('switch')
    SUB.append(sub)
    
    #repeat
    ns_al1=np.array([mean_ns_vertic_acc['1'],mean_ns_infreq_acc['1'],mean_ns_freq_acc['1']], np.float)
    ns_al1=np.nanmean(ns_al1)
    
    ACC.append(ns_al1)
    avACC.append(ns_al1)
    AL.append('AL1')
    SW.append('repeat')
    SUB.append(sub)
    ns_al2=np.array([mean_ns_vertic_acc['2'],mean_ns_infreq_acc['2'] ,mean_ns_freq_acc['2']],np.float)
    ns_al2=np.nanmean(ns_al2)
    ACC.append(ns_al2)
    avACC.append(ns_al2)
    AL.append('AL2')
    SW.append('repeat')
    SUB.append(sub)
    ns_al3=np.array([mean_ns_vertic_acc['3'], mean_ns_infreq_acc['3'] , mean_ns_freq_acc['3']], np.float)
    ns_al3=np.nanmean(ns_al3)
    
    ACC.append(ns_al3)
    avACC.append(ns_al3)
    AL.append('AL3')
    SW.append('repeat')
    SUB.append(sub)
     
    df=list(zip(SUB,ACC,AL,SW,avACC))
    
    df=pd.DataFrame(df)
    now = datetime.now()
    date = now.strftime("%m_%d_%Y")
    df.to_csv(datapath1 + 'mean_ACCs_long_allsubs_%s.csv'%(date), sep=',', index= False, header=True)

  
    table=[params,rts] 

    if sub ==1:
        with open(datapath1 + 'ACCs_allsubs_mg_beh_%s.csv' %(date) , 'w') as csvfile:
            writer = csv.writer(csvfile)
            [writer.writerow(r) for r in table]
    
    else:
        with open(datapath1 + 'ACCs_allsubs_mg_beh_%s.csv' %(date), 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(rts) 
   
    
    '''#####            SENSITIVITY ANALYSIS            #####'''


    '''sum up the rates of hits, misses etc.'''
    for key in SA_even:
        SA_even[key]=sum(SA_even[key]) 
    
    for key in SA_uneven:
        SA_uneven[key]=sum(SA_uneven[key]) 
    
    print (SA_even)
    print (SA_uneven)
    
    
    hits1=SA_even['Hit1']+SA_uneven['Hit1']
    false1=SA_even['False_Alarm1']+SA_uneven['False_Alarm1']
    missed1=SA_even['Missed_Hit1']+SA_uneven['Missed_Hit1']
    correctR1=SA_even['Correct_Reject1']+SA_uneven['Correct_Reject1']
    
            
    hits2=SA_even['Hit2']+SA_uneven['Hit2']
    false2=SA_even['False_Alarm2']+SA_uneven['False_Alarm2']
    missed2=SA_even['Missed_Hit2']+SA_uneven['Missed_Hit2']
    correctR2=SA_even['Correct_Reject2']+SA_uneven['Correct_Reject2']
    
            
    hits3=SA_even['Hit3']+SA_uneven['Hit3']
    false3=SA_even['False_Alarm3']+SA_uneven['False_Alarm3']
    missed3=SA_even['Missed_Hit3']+SA_uneven['Missed_Hit3']
    correctR3=SA_even['Correct_Reject3']+SA_uneven['Correct_Reject3']
    
            
    hits4=SA_even['Hit4']+SA_uneven['Hit4']
    false4=SA_even['False_Alarm4']+SA_uneven['False_Alarm4']
    missed4=SA_even['Missed_Hit4']+SA_uneven['Missed_Hit4']
    correctR4=SA_even['Correct_Reject4']+SA_uneven['Correct_Reject4']
    

    '''use function dPrime to calculate sensitivity measures'''

    d_1_uneven=dPrime(SA_uneven['Hit1'],SA_uneven['Missed_Hit1'],SA_uneven['False_Alarm1'],SA_uneven['Correct_Reject1'])
    d_2_uneven=dPrime(SA_uneven['Hit2'],SA_uneven['Missed_Hit2'],SA_uneven['False_Alarm2'],SA_uneven['Correct_Reject2'])
    d_3_uneven=dPrime(SA_uneven['Hit3'],SA_uneven['Missed_Hit3'],SA_uneven['False_Alarm3'],SA_uneven['Correct_Reject3'])
    d_4_uneven=dPrime(SA_uneven['Hit4'],SA_uneven['Missed_Hit4'],SA_uneven['False_Alarm4'],SA_uneven['Correct_Reject4'])

    '''for even blocks'''

    d_1_even=dPrime(SA_even['Hit1'],SA_even['Missed_Hit1'],SA_even['False_Alarm1'],SA_even['Correct_Reject1'])
    d_2_even=dPrime(SA_even['Hit2'],SA_even['Missed_Hit2'],SA_even['False_Alarm2'],SA_even['Correct_Reject2'])
    d_3_even=dPrime(SA_even['Hit3'],SA_even['Missed_Hit3'],SA_even['False_Alarm3'],SA_even['Correct_Reject3'])
    d_4_even=dPrime(SA_even['Hit4'],SA_even['Missed_Hit4'],SA_even['False_Alarm4'],SA_even['Correct_Reject4'])
    
    '''for all blocks'''
    d_1=dPrime(hits1,missed1,false1,correctR1)
    d_2=dPrime(hits2,missed2,false2,correctR2)
    d_3=dPrime(hits3,missed3,false3,correctR3)
    d_4=dPrime(hits4,missed4,false4,correctR4)
    
    
    c1=np.nanmean([d_1_even['c'],d_1_uneven['c']])
    c2=np.nanmean([d_2_even['c'],d_2_uneven['c']])
    c3=np.nanmean([d_3_even['c'],d_3_uneven['c']])
    c4=np.nanmean([d_4_even['c'],d_4_uneven['c']])

    a1=np.nanmean([d_1_even['A'],d_1_uneven['A']])
    a2=np.nanmean([d_2_even['A'],d_2_uneven['A']])
    a3=np.nanmean([d_3_even['A'],d_3_uneven['A']])
    a4=np.nanmean([d_4_even['A'],d_4_uneven['A']])
    
    pc1=np.nanmean([d_1_even['Pc'],d_1_uneven['Pc']])
    pc2=np.nanmean([d_2_even['Pc'],d_2_uneven['Pc']])
    pc3=np.nanmean([d_3_even['Pc'],d_3_uneven['Pc']])
    pc4=np.nanmean([d_4_even['Pc'],d_4_uneven['Pc']])
    
    all_d=['sub','c1','A1','c2','A2','c3','A3','c4','A4','Pc1','Pc2','Pc3','Pc4'
            ]
      
    vals=[sub,c1,a1,c2,a2,c3,a3,c4,a4,pc1,pc2,pc3,pc4]    
    
    table=[all_d,vals]
    
    if sub ==1:

        with open(datapath1 + 'averaged_SDTM_mg_behstudy_dec20_model_%s.csv'%(date), 'w') as csvfile:#write
            writer = csv.writer(csvfile)
            [writer.writerow(r) for r in table]
    else:
        with open(datapath1 + 'averaged_SDTM_mg_behstudy_dec20_model_%s.csv'%(date), 'a') as csvfile:#Ad
            writer = csv.writer(csvfile)
            writer.writerow(table[1]) 


    
 