# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 03:34:58 2020

@author: eleve2
"""
import math
import pyreadstat
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import readsav
from scipy.integrate import simps
import numpy as np
from scipy.signal import savgol_filter


result = readsav('158348data.sav')
time=result["time"]
dens=result["sumdens"]
fneg=result["fneg"]

dt=100 

nbins=100 

Tsize=len(time) 

tt=time[0:(Tsize-1):dt]

slidesize=int(Tsize/dt)

Esize=int((Tsize-slidesize)/dt)

E=np.linspace(0,1,Esize)

data=dens[3]

data=data[0:slidesize]


pp=np.histogram(data) 

pp=plt.hist(data,nbins,density=True)

pdf=[]

xx=[]

xx.append(pp[1])

pdf.append(savgol_filter(pp[0], 51, 3))

for k in range(0,Esize-2):
    
    temp1=data[k*dt:k*dt+slidesize]
    
    temp2=data[(k+1)*dt:(k+1)*dt+slidesize]
    #m1=max(temp1.max(),temp2.max())

    #m2=min(temp1.max(),temp2.max())

    p1=plt.hist(temp1,nbins,density=True)
    
    psm1=savgol_filter(p1[0], 51, 3)

    p2=plt.hist(temp2,nbins,density=True)
    
    psm2=savgol_filter(p2[0], 51, 3)

    integ=np.linspace(0,len(psm1),len(p1[1]))
    
    pdf.append(psm1)
    
    xx.append(p1[1])

    for l in range(0,len(integ)-2):
        if (psm1[l] >= 0) and (psm2[l] >= 0):
            integ[l]=psm1[l]*(math.log(psm2[l])-math.log(psm1[l]))**2/(dt**2) 

        else:
            integ[l]=0.000 #using the approximation 0*log(0)=0
            
    len(integ)
    len(xx[k])
    E[k]=simps(integ,xx[k])

#L=np.linspace(0,Esize)
"""
for n in range(1,(Esize-1)):
    L(n)=simps(T(0:n*dt:dt),double(sqrt(E(0:n))))
"""