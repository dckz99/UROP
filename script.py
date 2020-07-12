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
from scipy.integrate import trapz
import numpy as np
from scipy.signal import savgol_filter


result = readsav('158348data.sav')
time=result["time"]
dens=result["sumdens"]
fneg=result["fneg"]


m=10 #visualize m number of pdfs 
dataset=3 #there might be multiple datasets in the dens array
Tsize=len(time) #size of the dataset
res=10 #resolution
dt=100 #time discretization 
windowsize=int(Tsize/res) #size of the PDF window
nbins=100 #number of bins to use for the histogram affects accuracy of the PDF

tt=time[0:Tsize:dt] #time with dt increments

E=[] #initializing the array of the rate of change of information

L=[] #initializing the array of the information length function 

data=dens[dataset]

pp=np.histogram(data[0:windowsize],nbins,density=True) #first pdf

pdf=[]

xx=[]

xx.append(pp[1][:-1]) # x values for the pdf 

pdf.append(savgol_filter(pp[0], 51, 3)) # smoothening the pdf using a filter

k=0

while ((k+1)*dt+windowsize)<=Tsize:
    
    temp1=data[k*dt:k*dt+windowsize]
    
    temp2=data[(k+1)*dt:(k+1)*dt+windowsize]

    p1=np.histogram(temp1,nbins,density=True)
    
    psm1=savgol_filter(p1[0], 51, 3)

    p2=np.histogram(temp2,nbins,density=True)
    
    psm2=savgol_filter(p2[0], 51, 3)

    integ=[]
    
    pdf.append(psm2)
    
    xx.append(p2[1][:-1])

    for l in range(0,len(psm2)):
        if (psm1[l] > 0) and (psm2[l] > 0):
            integ.append(psm1[l]*(math.log(psm2[l])-math.log(psm1[l]))**2/(dt**2)) 

        else:
            integ.append(0.000) #using the approximation 0*log(0)~0

            
    E.append(math.sqrt(simps(integ,xx[k])))
    
    k=k+1

n=0

while ((n+1)*dt<=Tsize-windowsize):
    L.append(trapz(E[0:n+1],tt[0:(n+1)]))
    n=n+1

for i in range(0,m):
    plt.plot(xx[int(len(pdf)/m)*i+1],pdf[int(len(pdf)/m)*i+1],label='PDF at time '+str(int(len(pdf)/m)*i+1))
    plt.legend(loc='upper right')

    