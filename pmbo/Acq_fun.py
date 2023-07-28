#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import torch
from scipy.stats import norm


# In[2]:


def ei(y_pred, ybest):
    with torch.no_grad():
        mean = y_pred.mean.numpy()
        std = y_pred.stddev.numpy()
        y_max = ybest
        z = (mean-y_max) / std
        out = (mean-y_max) * norm.cdf(z) + std * norm.pdf(z)
    return out
def ucb(y_pred, y_train):
    out=[]
    with torch.no_grad():
        mean = y_pred.mean.numpy()
        std = y_pred.stddev.numpy()
        n_sample = y_train.numpy().shape[0]
        k=np.sqrt(np.log(n_sample)/n_sample)
        for i in range (0,n_sample):
            out_ = mean[i]+k* std[i]
            out.append (out_)
    return np.array(out)

