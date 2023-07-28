#!/usr/bin/env python
# coding: utf-8

# In[20]:


import math
import numpy as np
import torch
from copy import copy
import pandas as pd


# In[21]:


def centering(X):
    stdX = np.std(X, 0)
    index = np.where(stdX != 0)
    return (X[:, index[0]] - np.mean(X[:, index[0]], 0)) / stdX[index[0]]


# In[22]:


def cover_torch (data_array):
    if type(data_array)==torch.Tensor:
        temp_array=data_array
    else:
        temp_array=torch.tensor(data_array)        
    return temp_array
def cover_numpy (data_array):
    if type(data_array)==torch.Tensor:
        temp_array=data_array.numpy()
    else:
        temp_array=data_array
    return temp_array


# In[29]:


def Initialdata_Preparation (X,Xinfo,y,ntrain,natom_layer,n_job,seed):
    Xdata=np.array(X)
    Ydata = np.array(y)
    ndata = len(Ydata)
    nremain = ndata - ntrain
    np.random.seed(seed)
    dataset = np.random.permutation(ndata)
    a1data = np.empty(ntrain, dtype=int)
    a2data = np.empty(nremain, dtype=int)
    a1data[:] = dataset[0:ntrain]
    a2data[:] = dataset[ntrain:ndata]
    # info for the initial training set
    Xtrain = np.ndarray(shape=(ntrain, natom_layer), dtype=float)
    Xtraininfo = np.chararray(ntrain, itemsize=100)
    Ytrain = np.empty(shape=(ntrain, n_job), dtype=float)
    Xtrain[0:ntrain,:] = Xdata[a1data,:]
    Xtraininfo[0:ntrain] = Xinfo[a1data]
    Ytrain[0:ntrain,:] = Ydata[a1data,:]    
    #Ytrain[:,0:ntrain] = Ydata[:,a1data]  
    Xremain= np.ndarray(shape=(nremain, natom_layer), dtype=float)
    Xremaininfo = np.chararray(nremain, itemsize=100)
    Yremain = np.ndarray(shape=(nremain,n_job), dtype=float)
    Xremain[0:nremain, :] = Xdata[a2data, :]
    Xremaininfo[0:nremain] = Xinfo[a2data]
    Yremain[0:nremain,:] = Ydata[a2data,:]
    XtraininfoIni = Xtraininfo
    XtraininfoIni=np.array([x.decode() for x in XtraininfoIni])
    XremaininfoIni = Xremaininfo
    XremaininfoIni=np.array([x.decode() for x in XremaininfoIni])    
    return Xtrain,XtraininfoIni,Ytrain,Xremain,XremaininfoIni,Yremain


# In[24]:


def Search_polymer (polymer,polymerID,Threshold):
    wellpolymer=[]
    wellpolymerid=[]
    for i in range(len(polymer)):
        if (polymer[i][0]>=Threshold [0] and polymer[i][1]>=Threshold [1] and polymer[i][2]>=Threshold [2]):
            wellpolymer.append (polymer[i])
            wellpolymerid.append (polymerID[i])
    wellpolymer=np.array(wellpolymer)
    wellpolymerid=np.array(wellpolymerid)
    return wellpolymer,wellpolymerid   


# In[25]:


def matching_smiles(pdata,PID,IDNAME="monomer_ID",SMILENAME="smiles"):
    Datafarm= np.vstack((pdata[IDNAME].values,df[SMILENAME].values)).T
    for i in range(len(Datafarm)):
        if Datafarm[i][0]==PID:
            smi=Datafarm[i][1]
    return smi


# In[26]:


def Simulation_direct (ydata,location):
    yopt=ydata[location]
    return yopt


# In[27]:


def paretoSearch(capP,search='min'):
    # Non-dominated sorting
    paretoIdx=[]
    F0 = []
    for i,p in enumerate(capP):
        Sp = []
        nps = 0
        for j,q in enumerate(capP):
            if i!=j:
                if search=='min':
                    compare = p < q
                elif search=='max':
                    compare = p > q
                if any(compare):
                    Sp.append(q)
                else: 
                    nps+=1
        if nps==0:
            paretoIdx.append(i)
            prank = 1
            F0.append(p.tolist())
    F0 = np.array(F0)
    return F0, paretoIdx
def paretoOpt(capP, metric='crowdingDistance',opt='min'):
    if capP.shape[0]<=1000:
        F0, paretoIdx = paretoSearch(capP, search=opt)
    else:
        n_parts = int(capP.shape[0]//1000.)
        rem = capP.shape[0] % 1000.  
        FList = [] 
        paretoIdxList = []
        for i in range(n_parts):
            Fi, paretoIdxi = paretoSearch(capP[1000*i:1000*(i+1)], search=opt)
            FList.append(Fi)
            ar_paretoIdxi = np.array(paretoIdxi)+1000*i
            paretoIdxList.append(ar_paretoIdxi.tolist())  
        if rem>0:
            Fi, paretoIdxi = paretoSearch(capP[1000*n_parts-1:-1], search=opt)
            FList.append(Fi)
            ar_paretoIdxi = np.array(paretoIdxi)+1000*n_parts
            paretoIdxList.append(ar_paretoIdxi.tolist())  
            
        F1 = np.concatenate(FList)
        paretoIdx1=np.concatenate(paretoIdxList)
        F0, paretoIdxTemp = paretoSearch(F1, search=opt)
        paretoIdx=[]
        for a in paretoIdxTemp:
            matchingArr = np.where(capP==F1[a])[0]
            counts = np.bincount(matchingArr)
            pt = np.argmax(counts)
            paretoIdx.append(pt)
   
    m=F0.shape[-1]
    l = len(F0)
    ods = np.zeros(np.max(paretoIdx)+1)
    if metric == 'crowdingDistance':
        infi = 1E6
        for i in range(m):
            order = []
            sortedF0 = sorted(F0, key=lambda x: x[i])
            for a in sortedF0: 
                matchingArr = np.where(capP==a)[0]
                counts = np.bincount(matchingArr)
                o = np.argmax(counts)
                order.append(o)
            ods[order[0]]=infi
            ods[order[-1]]=infi
            fmin = sortedF0[0][i]
            fmax = sortedF0[-1][i]
            for j in range(1,l-1):
                ods[order[j]]+=(capP[order[j+1]][i]-capP[order[j-1]][i])/(fmax-fmin)
        # Impose criteria on selecting pareto points
        if min(ods[np.nonzero(ods)])>=infi:
            bestIdx = np.argmax(ods)
        else:
            if l>2: # if there are more than 2 pareto points, pick inner points with largest crowding distance (i.e most isolated)
                tempOds=copy(ods)
                for i,a in enumerate(tempOds):
                    if a>=infi: tempOds[i]=0.
                bestIdx = np.argmax(tempOds)
            else: #pick pareto point with lower index
                bestIdx = np.argmax(ods)
    elif metric == 'euclideanDistance':  # To the hypothetical point of the current data
        for i in range(m):
            order = []
            sortedF0 = sorted(F0, key=lambda x: x[i])
            for a in sortedF0:
                matchingArr = np.where(capP==a)[0]
                counts = np.bincount(matchingArr)
                o = np.argmax(counts)
                order.append(o)          
            fmin = sortedF0[0][i]
            fmax = sortedF0[-1][i]
            for j in range(0,l):
                ods[order[j]]+=((capP[order[j]][i]-fmax)/(fmax-fmin))**2
        ods = np.sqrt(ods)
        for i,a in enumerate(ods):
            if a!=0: print(i,a)
        bestIdx = np.where(ods==np.min(ods[np.nonzero(ods)]))[0][0]
    return paretoIdx,bestIdx    

