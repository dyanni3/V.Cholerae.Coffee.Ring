#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 10:48:34 2016
Bowl Shape
@author: dyanni3
"""

import numpy as np
from matplotlib import pyplot as plt

#%% Initialization -------------
x=(50*np.zeros(100))#+np.random.random(100) #cross section of agar surface
x[25]=1; x[75]=1

#%% Update rules --------------
def update(x):
    #s=np.random.choice(np.arange(1,99))
    a=x.copy()
    for s in range(25,76):
        dyr=a[s+1]-a[s]; dyl=a[s-1]-a[s]
        pitchr=max(.5-.25*dyr,0); pitchl=max(.5-.25*dyl,0)
        x[s+1]+=pitchr
        x[s-1]+=pitchl
        x[s]-=pitchr+pitchl
    a=x.copy()
    for s in range(1,len(x)-1):
        x[s-1]+=.2*(a[s]-a[s-1]); 
        x[s+1]+=.2*(a[s]-a[s+1]);
        x[s]-=.2*(a[s]-a[s+1])
        x[s]-=.2*(a[s]-a[s-1])
    return
#%% plotting ---------
plt.plot(x)
#%% Update rules for exponential growth style
def xUpdate(x):
    for s in range(1,99):
        x[s]+=.001*x[s] #exponential growth
        if np.random.random()>.1: #avalanche
            x[s+1],x[s-1]=x[s+1]+.5*(x[s]-x[s+1]),x[s-1]+.5*(x[s]-x[s-1])