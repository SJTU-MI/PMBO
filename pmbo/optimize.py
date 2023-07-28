#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from copy import copy


# In[1]:


from utility import cover_torch,cover_numpy,paretoSearch,paretoOpt,Simulation_direct
from gp import gp_model,gp_predict
from Acq_fun import ei


# In[2]:


def MBO_opt(Xtrain,XtraininfoIni,Ytrain,Xremain,XremaininfoIni,Yremain,nobj,n_layer,Multiprobe,X_layer,ybest,metr):
    modelList = []
    likehoodList = []
    yt_predList = []
    TXtrain=cover_torch (Xtrain)
    TYtrain=cover_torch (Ytrain)
    TXremain=cover_torch (Xremain)
    TYremain=cover_torch (Yremain)
    for i in range(nobj):
        model,likelihood=gp_model (TXtrain[:,X_layer[0][i]:(X_layer[0][i]+X_layer[0][i+1])],TYtrain[:,i],5000,error=1e-10)
        modelList.append(model)
        likehoodList.append(likehoodList)
        yt_pred=gp_predict(model,likelihood,TXremain[:,X_layer[0][i]:(X_layer[0][i]+X_layer[0][i+1])]) 
        yt_predList.append(yt_pred)
        if i==0:
            yt_mean=yt_predList[i].mean.numpy()
            EI=ei(yt_pred, ybest[i])
        else:
            yt_mean=np.vstack((yt_mean,yt_predList[i].mean.numpy()))
            EI=np.vstack((EI,ei(yt_pred, ybest[i])))
        del model,likelihood
    yt_mean=yt_mean.T
    EI=EI.T
    EInew=copy(EI)
    EIrank=np.arange(len(EInew))
    RECloc=[]
    for ii in range (Multiprobe):
        if ii==0:
            _, EIbestloc = paretoOpt(EInew,metric=metr,opt='max') 
            RECloc.append(EIrank[EIbestloc])
        else:
            EInew=np.delete(EInew, EIbestloc, 0)
            EIrank=np.delete(EIrank, EIbestloc, 0)
            _, EIbestloc = paretoOpt(EInew,metric=metr,opt='max') 
            RECloc.append(EIrank[EIbestloc])
    Xexplored=Xremain[RECloc]
    Yexplored=Simulation_direct(Yremain,RECloc)
    Xexploredinfo=XremaininfoIni[RECloc]
    Xtrainnew=np.append(Xtrain,Xremain[RECloc]).reshape(-1,n_layer)
    XtraininfoIninew=np.append(XtraininfoIni,XremaininfoIni[RECloc])
    Ytrainnew=np.append(Ytrain, Yexplored).reshape(-1,nobj)
    Xremainnew= np.delete(Xremain, RECloc, 0)
    XremaininfoIninew=np.delete(XremaininfoIni,RECloc, 0)
    Yremainnew=np.delete(Yremain,RECloc, 0)
    _, ybestloc = paretoOpt(Ytrainnew,metric=metr,opt='max') 
    ybestnew=Ytrainnew[ybestloc]
    ybestidnew=XtraininfoIninew[ybestloc]
    result=XremaininfoIni[RECloc]
    del EI,EInew,EIrank,RECloc
    return  Xtrainnew,XtraininfoIninew,Ytrainnew,Xremainnew,XremaininfoIninew,Yremainnew,Xexplored,Xexploredinfo,Yexplored,ybestnew,ybestidnew,result


# In[ ]:




