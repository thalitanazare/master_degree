#%%
# Artigo Aguirre1995: Retrieving Dynamical Invariants from Chaotic Data Using NARMAX Models

"""
Usando o somatório de Kahan
A referência é o modelo da literatura simulado com maior precisão
usando o mpmath. Todas variáveis devem ser decladas com esse pacote, assim 
como as operações matemáticas
"""
#==============================================================================
#=======================...BIBLIOTECAS...======================================
#==============================================================================
import numpy as np
from numpy import random
import pandas as pd
import matplotlib.pyplot as plt
import math as mt
from decimal import *
plt.rcParams["font.family"] = "Times New Roman" 
import mpmath as mp
from mpmath import *
from sympy import Symbol, solveset, S
from sklearn.linear_model import LinearRegression
import random as r
import scipy.stats as stats
#==============================================================================
#=======================......KAHAN......======================================
#==============================================================================

def kahansum(input):
    summ = c = 0
    for num in input:
        y = num - c
        t = summ + y
        c = (t - summ) - y
        summ = t
    return summ


#==============================================================================
#=======================.....ENTRADAS....======================================
#==============================================================================
mp.dps=34





#------------------------------------------------------------------------------
# DADOS MATLAB
dados_normal=np.loadtxt("xmackeyglass_normal.txt")
dados_kahan=np.loadtxt("xmackeyglass_kahan.txt")
xmg_normal=dados_normal
xmg_kahan=dados_kahan
#-------------------------
N=5000
mp.dps=34
#condição inicial usando mpf
#x=[mp.mpf(x0[0]),mp.mpf(x0[1]),mp.mpf(x0[2]),mp.mpf(x0[3]),mp.mpf(x0[4]),mp.mpf(x0[5]),mp.mpf(x0[6]),mp.mpf(x0[7]),mp.mpf(x0[8]),mp.mpf(x0[9])]
#x0=[0.32185972,0.40354956,0.16261371,0.19585123,0.83952246,0.20231749,0.15468285,0.50401686,0.42225821,0.34911833]
#x0=[0.3009654801388982, 0.3004361618666274, 0.30062664829086677, 0.30030102619842547, 0.300507242983829, 0.3003858662588449, 0.30035091048877016, 0.30058507410740537, 0.3005842517929702, 0.30090420177084776]
#xref=[mp.mpf(x0[0]),mp.mpf(x0[1]),mp.mpf(x0[2]),mp.mpf(x0[3]),mp.mpf(x0[4]),mp.mpf(x0[5]),mp.mpf(x0[6]),mp.mpf(x0[7]),mp.mpf(x0[8]),mp.mpf(x0[9])]
h=10
cond=[]
r.seed(2)
#condição inicial com maior precisãão
for i in range(1,11): 
    cond=np.append(cond, r.uniform(0.01,0.01))

x0=[cond[0],cond[1],cond[2],cond[3],cond[4],cond[5],cond[6],cond[7],cond[8],cond[9]]
xref = [mp.mpf(x0[0]),mp.mpf(x0[1]),mp.mpf(x0[2]),mp.mpf(x0[3]),mp.mpf(x0[4]),mp.mpf(x0[5]),mp.mpf(x0[6]),mp.mpf(x0[7]),mp.mpf(x0[8]),mp.mpf(x0[9])]
# Modelo do sistema referência


for k in range(h,N):
    xref = np.append(xref, mp.mpf(0.24662*10) *xref[k-1] - mp.mpf(0.16423*10)*xref[k-2]\
                    + mp.mpf(0.60992)*xref[k-3] + mp.mpf(0.073012)*xref[k-5]**2*xref[k-10]**2\
                    + mp.mpf(0.38566)*xref[k-3]*xref[k-10]\
                    + mp.mpf(0.66999)*xref[k-1]*xref[k-10]**2\
                    + mp.mpf(0.88364)*xref[k-1]**3\
                    - mp.mpf(0.67300)*xref[k-4]*xref[k-10]**2\
                    - mp.mpf(0.11929*10)*xref[k-1]**2\
                    - mp.mpf(0.050451)*xref[k-4]*xref[k-5] - mp.mpf(0.24765)*xref[k-1]**4\
                    + mp.mpf(0.42081)*xref[k-4]*xref[k-9]*xref[k-10]**2\
                    - mp.mpf(0.70406)*xref[k-1]*xref[k-10]**3\
                    - mp.mpf(0.14089)*xref[k-5]*xref[k-8]**2\
                    + mp.mpf(0.14807)*xref[k-1]*xref[k-7]*xref[k-10]   )
#------------------------------------------------------------------------------

# Modelo do sistema tradicional
xnumpy=x0

for k in range(h,N):
    xnumpy = np.append(xnumpy, 0.24662*10*xnumpy[k-1] - 0.16423*10*xnumpy[k-2]\
                    + 0.60992*xnumpy[k-3] + 0.073012*xnumpy[k-5]**2*xnumpy[k-10]**2\
                    + 0.38566*xnumpy[k-3]*xnumpy[k-10]\
                    + 0.66999*xnumpy[k-1]*xnumpy[k-10]**2\
                    + 0.88364*xnumpy[k-1]**3\
                    - 0.67300*xnumpy[k-4]*xnumpy[k-10]**2\
                    - 0.11929*10*xnumpy[k-1]**2\
                    - 0.050451*xnumpy[k-4]*xnumpy[k-5] - 0.24765*xnumpy[k-1]**4\
                    + 0.42081*xnumpy[k-4]*xnumpy[k-9]*xnumpy[k-10]**2\
                    - 0.70406*xnumpy[k-1]*xnumpy[k-10]**3\
                    - 0.14089*xnumpy[k-5]*xnumpy[k-8]**2\
                    + 0.14807*xnumpy[k-1]*xnumpy[k-7]*xnumpy[k-10]   )      
    
xkahan=x0
for k in range(h,N):
    xkahan = np.append(xkahan, kahansum([0.24662*10*xkahan[k-1],- 0.16423*10*xkahan[k-2],\
                     0.60992*xkahan[k-3], + 0.073012*xkahan[k-5]**2*xkahan[k-10]**2,\
                     0.38566*xkahan[k-3]*xkahan[k-10],\
                     0.66999*xkahan[k-1]*xkahan[k-10]**2,\
                     0.88364*xkahan[k-1]**3,\
                    - 0.67300*xkahan[k-4]*xkahan[k-10]**2,\
                    - 0.11929*10*xkahan[k-1]**2,\
                    - 0.050451*xkahan[k-4]*xkahan[k-5], - 0.24765*xkahan[k-1]**4,\
                     0.42081*xkahan[k-4]*xkahan[k-9]*xkahan[k-10]**2,\
                    - 0.70406*xkahan[k-1]*xkahan[k-10]**3,\
                    - 0.14089*xkahan[k-5]*xkahan[k-8]**2,\
                     0.14807*xkahan[k-1]*xkahan[k-7]*xkahan[k-10]])  ) 


#==============================================================================
#===================== EXPOENTE DE LYAPUNOV ===================================
#---------------------------- + NUMPY +----------------------------------------

# PARÂMETROS PARA PONTOS DA LINEARREGRESSION E LBE
start=50 #ponto que começa calcular o LBE (de preferência o mesmo de h)
vec1=150 #primeiro ponto da reta
vec2=1490 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEpy=abs(xref[start:N]-xnumpy[start:N])/2 #calcula o LBE
LBEpy=np.array(LBEpy,dtype=np.float) #transformando para array tipo np
LBEnumpy=np.log2(LBEpy) #gera o LBE no gráfico de log2

LLEnumpy=LinearRegression()
LLEnumpy.fit(vecN.reshape(-1,1),LBEnumpy[vec1:vec2])
LLEpy_coef=LLEnumpy.coef_
LLEpy_inter=LLEnumpy.intercept_

Ncnumpy=int(abs(np.log10(1.388/2)+np.log2(10**16))/LLEpy_coef)

print("-----------------------")
print("LLE numpy")
print(LLEpy_coef)
print("Tempo numpy")
print(Ncnumpy)


# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEnumpy,color='black',linewidth=0.75)
plt.plot(vecN,(vecN*LLEpy_coef+LLEpy_inter),color='red',linewidth=3)
plt.axvline(x=Ncnumpy, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEnumpy)-3, color='black',linestyle="solid")
plt.text(Ncnumpy+1, -53, 'N$_c$= %s'%Ncnumpy,fontsize=16, color="red" )
plt.text(2000, -15, '%0.4fn'%LLEpy_coef,fontsize=17 )
plt.text(2670, -15, '%0.2f'%LLEpy_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEnumpy)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{numpy}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('mg_lbe_numpy.pdf', format='pdf')

#---------------------------- + NUMPY KAHAN +-----------------------------------
start=10
vec1=50 #primeiro ponto da reta
vec2=1800#último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEkhpy=abs(xref[start:N]-xkahan[start:N])/2 #calcula o LBE
LBEkhpy=np.array(LBEkhpy,dtype=np.float) #transformando para array tipo np
LBEkahanpy=np.log2(LBEkhpy) #gera o LBE no gráfico de log2

LLEkahanpy=LinearRegression()
LLEkahanpy.fit(vecN.reshape(-1,1),LBEkahanpy[vec1:vec2])
LLEkhpy_coef=LLEkahanpy.coef_
LLEkhpy_inter=LLEkahanpy.intercept_


Nckahanpy=int(abs(np.log10(max(LBEkhpy))+np.log2(10**16))/LLEkahanpy.coef_)

print("-----------------------")
print("LLE numpy KAHAN")
print(LLEkhpy_coef)
print("Tempo numpy KAHAN")
print(Nckahanpy)

# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEkahanpy,color='black',linewidth=0.75)
plt.plot(vecN,(vecN*LLEkhpy_coef+LLEkhpy_inter),color='red',linewidth=3)
plt.axvline(x=Nckahanpy, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEkahanpy)-3, color='black',linestyle="solid")
plt.text(Nckahanpy+1, -63, 'N$_c$= %s'%Nckahanpy,fontsize=16, color="red" )
plt.text(2000, -19, '%0.4fn'%LLEkhpy_coef,fontsize=17 )
plt.text(2670, -19, '%0.2f'%LLEkhpy_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEkahanpy)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{kahan}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('mg_lbe_kahanpy.pdf', format='pdf')
#==============================================================================
#CÁLCULO NRMSE
#==============================================================================

n=int(1800/10)
#python
NRMSE1=[]
for i in range(1,11):
    NRMSE1=np.append(NRMSE1, np.sqrt(np.sum((xref[0:i*n] - xnumpy[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xnumpy[0:i*n])**2)))))

NRMSE2=[]
for i in range(1,11):
    NRMSE2=np.append(NRMSE2, np.sqrt(np.sum((xref[0:i*n] - xkahan[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xkahan[0:i*n])**2)))))

print('============================================================')
print('PYTHON POR PARTES-----> NRMSE REFERÊNCIA X NUMPY')
NRMSE_numpy=np.array(NRMSE1.tolist(),dtype=np.float64) #transformar mpf to numpy array
NRMSE_kahanpy=np.array(NRMSE2.tolist(),dtype=np.float64) #transformar mpf to numpy array
title=[" ","NUMPY","PYTHON KAHAN"]
it=np.linspace(n,1800, num=10)
it=it.astype(int)
#NRMSE=[[NRMSE_numpy] [NRMSE_kahanpy]]
NRMSE_2=pd.DataFrame({ 'numpy': NRMSE_numpy,
                    'python kahan': NRMSE_kahanpy})
NRMSE_2.style
NRMSE_2=NRMSE_2.set_index([pd.Index(it)])
#NRMSE.style
print(NRMSE_2)

plt.figure(figsize=(8, 6))
plt.plot(xref,color='black',linewidth=2, alpha=0.8,marker = None)
plt.plot(xnumpy,color='red',linewidth=1,marker=None)
plt.axvline(x=Ncnumpy, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Ncnumpy+3,0.05, 'N$_c$= %s'%Ncnumpy,fontsize=16 )
plt.axvline(x=1500, color='black',linestyle="solid")
plt.axhline(y=0, color='black',linestyle="solid")
plt.axis([1500, 2100, 0, 1.5])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{numpy}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('mg_ref_numpy.pdf', format='pdf')

plt.figure(figsize=(8, 6))
plt.plot(xref,color='black',linewidth=2, alpha=0.8,marker = None)
plt.plot(xkahan,color='red',linewidth=1,marker=None)
plt.axvline(x=Nckahanpy, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Nckahanpy+3,0.05, 'N$_c$= %s'%Nckahanpy,fontsize=16 )
plt.axvline(x=1500, color='black',linestyle="solid")
plt.axhline(y=0, color='black',linestyle="solid")
plt.axis([1500, 2100, 0, 1.5])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N kahan',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{kahan}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('mg_ref_kahanpy.pdf', format='pdf')


#%%
#==============================================================================
#===================== EXPOENTE DE LYAPUNOV ===================================
#---------------------------- + MATALB +---------------------------------------


# PARÂMETROS PARA PONTOS DA LINEARREGRESSION E LBE
start=100 #ponto que começa calcular o LBE (de preferência o mesmo de h)
vec1=100 #primeiro ponto da reta
vec2=1450 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEmat=abs(xref[start:N]-xmg_normal[start:N])/2 #calcula o LBE
LBEmat=np.array(LBEmat,dtype=np.float64) #transformando para array tipo np
LBEmatlab=np.log2(LBEmat) #gera o LBE no gráfico de log2

LLEmatlab=LinearRegression()
LLEmatlab.fit(vecN.reshape(-1,1),LBEmatlab[vec1:vec2])
LLEmat_coef=LLEmatlab.coef_
LLEmat_inter=LLEmatlab.intercept_

Ncmatlab=int(abs(np.log10(1.3887/2)+np.log2(10**16))/LLEmat_coef)

print("-----------------------")
print("LLE Matlab")
print(LLEmat_coef)
print("Tempo Matlab")
print(Ncmatlab)

# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEmatlab,color='black',linewidth=0.75)
plt.plot(vecN,(vecN*LLEmat_coef+LLEpy_inter),color='red',linewidth=3)
plt.axvline(x=Ncmatlab, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEmatlab)-3, color='black',linestyle="solid")
plt.text(Ncmatlab-4, -54, '%s'%Ncmatlab,fontsize=16, color="red" )
plt.text(2000, -15, '%0.4fn'%LLEmat_coef,fontsize=17 )
plt.text(2670, -15, '%0.2f'%LLEmat_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEmatlab)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{matlab}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('mg_lbe_matlab.pdf', format='pdf')
#---------------------------- + MATLAB KAHAN +-----------------------------------
star=10
vec1=0 #primeiro ponto da reta
vec2=1700 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEkhmat=abs(xref[start:N]-xmg_kahan[start:N])/2 #calcula o LBE
LBEkhmat=np.array(LBEkhmat,dtype=np.float64) #transformando para array tipo np
LBEkahanmat=np.log2(LBEkhmat) #gera o LBE no gráfico de log2

LLEkahanmat=LinearRegression()
LLEkahanmat.fit(vecN.reshape(-1,1),LBEkahanmat[vec1:vec2])
LLEkhmat_coef=LLEkahanmat.coef_
LLEkhmat_inter=LLEkahanmat.intercept_


Nckahanmat=int(abs(np.log10(1.119/2)+np.log2(10**16))/LLEkahanmat.coef_)

print("-----------------------")
print("LLE matlab KAHAN")
print(LLEkhmat_coef)
print("Tempo matlab KAHAN")
print(Nckahanmat)



# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEkahanmat,color='black',linewidth=0.75)
plt.plot(vecN,(vecN*LLEkhmat_coef+LLEkhmat_inter),color='red',linewidth=3)
plt.axvline(x=Nckahanmat, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEkahanmat)-3, color='black',linestyle="solid")
plt.text(Nckahanmat-4, -57, 'N$_c$= %s'%Nckahanmat,fontsize=16, color="red" )
plt.text(2000, -17, '%0.4fn'%LLEkhmat_coef,fontsize=17 )
plt.text(2670, -17, '%0.2f'%LLEkhmat_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEkahanmat)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{kahan}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('mg_lbe_kahanmat.pdf', format='pdf')


#==============================================================================
#=======================-----NRMSE-------======================================
#==============================================================================
n=int(1800/10)
#----------------------------matlab--------------------------------------------
#matlab
NRMSE_normal=[]
for i in range(2,12):
    NRMSE_normal=np.append(NRMSE_normal, np.sqrt(np.sum((xref[0:i*n] - xmg_normal[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xmg_normal[0:i*n])**2)))))


NRMSE_kahan=[]
for i in range(2,12):
    NRMSE_kahan=np.append(NRMSE_kahan, np.sqrt(np.sum((xref[0:i*n] - xmg_kahan[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xmg_kahan[0:i*n])**2)))))



print('============================================================')
print('MATALAB POR PARTES-----> NRMSE REFERÊNCIA X COM KAHAN')
NRMSE_kahan=np.array(NRMSE_kahan.tolist(),dtype=np.float64) #transformar mpf to numpy array
NRMSE_normal=np.array(NRMSE_normal.tolist(),dtype=np.float64) #transformar mpf to numpy array
title=[" ","MATLAB","MATLAB KAHAN"]
it=np.linspace(n,1800, num=10)
it=it.astype(int)
#NRMSE=[[NRMSE_numpy] [NRMSE_kahanpy]]
NRMSE_mat=pd.DataFrame({ 'matlab': NRMSE_normal,
                    'matlab kahan': NRMSE_kahan})
NRMSE_mat.style
NRMSE_mat=NRMSE_mat.set_index([pd.Index(it)])
#NRMSE.style
print(NRMSE_mat)


plt.figure(figsize=(8, 6))
plt.plot(xref,color='black',linewidth=2, alpha=0.8,marker = None)
plt.plot(xmg_normal,color='red',linewidth=1,marker=None)
plt.axvline(x=Ncmatlab, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Ncmatlab+3,0.05, 'N$_c$= %s'%Ncmatlab,fontsize=16 )
plt.axvline(x=1500, color='black',linestyle="solid")
plt.axhline(y=0, color='black',linestyle="solid")
plt.axis([1500, 2100, 0, 1.5])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{matlab}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('mg_ref_matlab.pdf', format='pdf')

plt.figure(figsize=(8, 6))
plt.plot(xref,color='black',linewidth=2, alpha=0.8,marker = None)
plt.plot(xmg_kahan,color='red',linewidth=1,marker=None)
plt.axvline(x=Nckahanmat, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Nckahanmat+3,0.05, 'N$_c$= %s'%Nckahanmat,fontsize=16 )
plt.axvline(x=1500, color='black',linestyle="solid")
plt.axhline(y=0, color='black',linestyle="solid")
plt.axis([1500, 2000, 0, 1.5])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{kahan}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('mg_ref_kahanmat.pdf', format='pdf')



