# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 15:46:42 2021

@author: thali
"""

import numpy as np #Usado para construir os dados de entrada do algoritmo
from numpy import random
import matplotlib.pyplot as plt #Gerar gráficos
import pandas as pd #Exibit tabelas organizadas com DataFrame
import math as mt
# Pacotes para intervalos
import mpmath as mp
from mpmath import *
from decimal import *
plt.rcParams["font.family"] = "Times New Roman" 

r=4 #parametro de bifurcação
x0=(mt.sqrt(3)+2)/4 #cond ini controlada
#x0= 0.4
x=[x0]
x1=x
N=100

precis=34
mp.dps=precis
rmp=mp.mpf('4') #parametro de bifurcação
x0mp=(mp.sqrt(3)+2)/4 #cond ini controlada
#x0= 0.4
xmp=[x0mp]
x1mp=xmp

for k in range (N):
    x = np.append(x, r*x[k]*(1-x[k])) #F

for k in range (N):
    xmp = np.append(xmp, rmp*xmp[k]*(1-xmp[k])) #Fmpf



plt.figure(figsize=(6,5))
plt.plot(x,color='black',linewidth=1.5,marker="x")
plt.plot(xmp,color='red',linewidth=2.5,alpha=1)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=-0.15, color='black',linestyle="solid")
plt.axis([0, 100, -0.15, 1.15])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=16)
plt.ylabel('F(x$_n$)',fontsize=16)
plt.grid(False)
plt.box(False)
plt.savefig('mapalog_mpf.pdf', format='pdf')

