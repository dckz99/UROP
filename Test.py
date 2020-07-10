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

Le ven. 10 juil. 2020 à 15:27, Adel Touzani <adeltouzani@gmail.com> a écrit :
IDL Routine slidescript with smoothing final:

T=time[value_locate(time,[0])+1:*] ;removing negative values for time
Y=Yval[value_locate(time,[0])+1:*] ;and the corresponding Y values
dt=100 ;dt discretisation of the time differential

nbindiv=100 ;divide nbins by this factor
smoothfact=10 ;smoothening factor

TSize=size(T,/n_elements) ;number of time points
tt=T(0:Tsize-1:dt) ;basically time but by increments of dt
YSize=size(Y,/n_elements) ;number of Y value points

slidesize=TSize/100
Esize=(Tsize-slidesize)/dt ;we need to make sure that we have a least a slidesize interval at the end of the algorithm
E=dblarr(Esize)

pp=double(histogram(Y(0:slidesize),nbins=slidesize/nbindiv,locations=xbinvals)) ;this will enable us to store the pdfs for later study
xx=xbinvals
pp=smooth(pp,size(pp,/n_elements)/smoothfact)
for k=0,Esize-2 do begin
k=long(k)

temp1=Y(k*dt:k*dt+slidesize) ;data to use for generating the first pdf
temp2=Y((k+1)*dt:(k+1)*dt+slidesize) ;data to use for the second
m1=max([max(temp1),max(temp2)])

m2=min([min(temp1),min(temp2)])

p1 = double(histogram(temp1,nbins=slidesize/nbindiv,locations=xbinvals,max=m1,min=m2))
psm1=smooth(p1,size(p1,/n_elements)/smoothfact)
p1t=int_tabulated(xbinvals,psm1)

psm1=psm1/p1t ;pdf at t


p2 = double(histogram(temp2,nbins=slidesize/nbindiv,locations=xbinvals,max=m1,min=m2))
psm2=smooth(p2,size(p2,/n_elements)/smoothfact)
p2t=int_tabulated(xbinvals,psm2)


psm2=psm2/p2t ;pdf at t+dt


int=dblarr(size(psm1,/n_elements)) ;the integrand
pp=[[pp],[psm1]]
xx=[[xx],[xbinvals]]
for l=0,(size(int,/n_elements)-2) do begin
if (psm1(l) ne 0) and (psm2(l) ne 0) then begin
int(l)=psm1(l)*(alog(psm2(l))-alog(psm1(l)))^2/(dt^2) ;
;the integrand is defined by the formula provided by Dr.Kim, it looks like some form of entropy to me which ;is to be ;expected from calcultaing the information length

endif else begin
int(l)=0.000 ;using the approximation 0*log(0)=0
endelse
endfor


E(k)=int_tabulated(xx(*,k),int)
endfor

L=dblarr(Esize)

for n=1,(Esize-1) do begin
n=long(n)
L(n)=int_tabulated(T(0:n*dt:dt),double(sqrt(E(0:n))))


endfor

end
