# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 15:46:42 2021

@author: thali
"""

import numpy as np #Usado para construir os dados de entrada do algoritmo
from numpy import random
import matplotlib.pyplot as plt #Gerar gráficos
import pandas as pd #Exibit tabelas organizadas com DataFrame
plt.rcParams["font.family"] = "Times New Roman" 

r=3.9 #parametro de bifurcação
x0=0.1 #cond ini controlada
#x0= 0.4
x=[x0]
y=[x0]
N=200



for k in range (N):
    x = np.append(x, r*x[k]*(1-x[k])) #F

for k in range (N):
    y = np.append(y, r*(y[k]*(1-y[k]))) #G

start=0
LBE=abs(y[start:N]-x[start:N])/2
graf_LBE=np.log10(LBE)

vec1=4
vec2=80
num1=vec2-vec1
vecN=np.linspace(vec1,vec2,num=num1)

LLE=np.polyfit(vecN,graf_LBE[vec1:vec2],1)
LLE1=LLE[0]
LLE2=LLE[1]

Ts=int(abs((np.log10(1/2)+16)/LLE[0]))

plt.figure(figsize=(8,4))
plt.plot(x,color='black',linewidth=1,marker="o")
plt.plot(y,color='red',linewidth=1,marker="x",alpha=0.8)
plt.axvline(x=40, color='black',linestyle="solid")
plt.axhline(y=-0.15, color='black',linestyle="solid")
plt.axis([40, 100, -0.15, 1.15])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=16)
plt.ylabel('F(x$_n$) , G(x$_n$)',fontsize=16)
plt.grid(False)
plt.box(False)
plt.savefig('mapalogEX.pdf', format='pdf')

#Numpy
plt.figure(figsize=(8, 6))
plt.plot(graf_LBE,color='black',linewidth=2,marker="o")
plt.plot(vecN,np.polyval(LLE,vecN),color='red',linewidth=4.3)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=-17, color='black',linestyle="solid")
plt.text(40, -16, '%0.5fn'%LLE1,fontsize=17 )
plt.text(65, -16, '%0.5f'%LLE2 ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, 160, -17, 1])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_{10}|\delta_{α,n}|$',fontsize=17)
plt.grid(False)
plt.box(False)
plt.savefig('mapalogLLE.pdf', format='pdf')
#%%
#Numpy
plt.figure(figsize=(8, 6))
plt.plot(graf_LBE,color='black',linewidth=2,marker="o")
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=-17, color='black',linestyle="solid")
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, 160, -17, 1])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_{10}|\delta_{α,n}|$',fontsize=17)
plt.grid(False)
plt.box(False)
plt.savefig('mapalogLBE.pdf', format='pdf')


#Numpy
plt.figure(figsize=(8, 6))
plt.plot(graf_LBE,color='black',linewidth=3,marker="o")
plt.plot(vecN,np.polyval(LLE,vecN),color='red',linewidth=4.3)
plt.axvline(x=Ts, color='red',linestyle='--',alpha=0.9,linewidth=1.5)
plt.text(Ts+2, -16, 'N$_c$ = %s'%Ts,fontsize=17 )
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=-17, color='black',linestyle="solid")
plt.text(25, -16, '%0.5fn'%LLE1,fontsize=17 )
plt.text(48, -16, '%0.5f'%LLE2 ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, 160, -17, 1])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_{10}|\delta_{α,n}|$',fontsize=17)
plt.grid(False)
plt.box(False)
plt.savefig('mapalogTS.pdf', format='pdf')