#%%
#rossler

"""
Usando o somatório de Kahan e usando horner
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
r.seed(10)
cond=[]
#condição inicial com maior precisãão
for i in range(1,6): #preciso de 5 cond iniciais e rodar 30 vezes
    cond=np.append(cond, r.uniform(0.1,0.10001))

x0=[cond[0],cond[1],cond[2],cond[3],cond[4]]
xref = [mp.mpf(x0[0]),mp.mpf(x0[1]),mp.mpf(x0[2]),mp.mpf(x0[3]),mp.mpf(x0[4])]
h=5
N=5000
# maior que isso dá NaN


# DADOS MATLAB
dados_normal=np.loadtxt("xrossler_normal.txt")
dados_kahan=np.loadtxt("xrossler_kahan.txt")
xrossler_normal=dados_normal
xrossler_kahan=dados_kahan


#==============================================================================
#=======================....SIMULAÇÃO....======================================
#==============================================================================

# Primeira equação sem modificação
for i in range(h,N):
    xref = np.append( xref, mp.mpf(0.1972*10)*xref[i-1] - mp.mpf(0.104*10)*xref[i-2] \
                  + mp.mpf(0.7456e-4)*xref[i-4]*xref[i-2]**3*xref[i-1] \
                  - mp.mpf(0.2053e-4)*xref[i-5]*xref[i-4]**4\
                  - mp.mpf(0.285e-4)*xref[i-5]*xref[i-1]**4\
                  + mp.mpf(0.2484e-4)*xref[i-3]**2*xref[i-2]**3\
                  + mp.mpf(0.1238e-2)*xref[i-2]*xref[i-1]**2\
                  + mp.mpf(0.4353e-4)*xref[i-5]**4\
                  + mp.mpf(0.2258e-2)*xref[i-5]*xref[i-2]*xref[i-1]**2\
                  + mp.mpf(0.3123e-4)*xref[i-4]**5\
                  + mp.mpf(0.7531e-3)*xref[i-1]**4\
                  - mp.mpf(0.2703e-2)*xref[i-3]**2*xref[i-1]**2\
                  - mp.mpf(0.7807e-3)*xref[i-1]**3\
                  - mp.mpf(0.7077e-4)*xref[i-3]**2*xref[i-2]**2*xref[i-1]\
                  - mp.mpf(0.3304e-3)*xref[i-3]*xref[i-2]**3\
                  - mp.mpf(0.8847e-2)*xref[i-5]*xref[i-1]\
                  + mp.mpf(0.7631e-2)*xref[i-4]*xref[i-1]\
                  - mp.mpf(0.3870e-4)*xref[i-5]**3*xref[i-1]**2\
                  + mp.mpf(0.4676e-3)*xref[i-3]**3*xref[i-1] )
    

xnumpy=x0
for i in range(h,N):
    xnumpy = np.append( xnumpy, 0.1972*10*xnumpy[i-1] - 0.104*10*xnumpy[i-2] \
                  + 0.7456e-4*xnumpy[i-4]*xnumpy[i-2]**3*xnumpy[i-1] \
                  - 0.2053e-4*xnumpy[i-5]*xnumpy[i-4]**4\
                  - 0.285e-4*xnumpy[i-5]*xnumpy[i-1]**4\
                  + 0.2484e-4*xnumpy[i-3]**2*xnumpy[i-2]**3\
                  + 0.1238e-2*xnumpy[i-2]*xnumpy[i-1]**2\
                  + 0.4353e-4*xnumpy[i-5]**4\
                  + 0.2258e-2*xnumpy[i-5]*xnumpy[i-2]*xnumpy[i-1]**2\
                  + 0.3123e-4*xnumpy[i-4]**5\
                  + 0.7531e-3*xnumpy[i-1]**4\
                  - 0.2703e-2*xnumpy[i-3]**2*xnumpy[i-1]**2\
                  - 0.7807e-3*xnumpy[i-1]**3\
                  - 0.7077e-4*xnumpy[i-3]**2*xnumpy[i-2]**2*xnumpy[i-1]\
                  - 0.3304e-3*xnumpy[i-3]*xnumpy[i-2]**3\
                  - 0.8847e-2*xnumpy[i-5]*xnumpy[i-1]\
                  + 0.7631e-2*xnumpy[i-4]*xnumpy[i-1]\
                  - 0.3870e-4*xnumpy[i-5]**3*xnumpy[i-1]**2\
                  + 0.4676e-3*xnumpy[i-3]**3*xnumpy[i-1] )
        
        
xkahan=x0
for i in range(h,N):
    xkahan = np.append( xkahan, kahansum([0.1972*10*xkahan[i-1], - 0.104*10*xkahan[i-2], \
                  + 0.7456e-4*xkahan[i-4]*xkahan[i-2]**3*xkahan[i-1], \
                  - 0.2053e-4*xkahan[i-5]*xkahan[i-4]**4,\
                  - 0.285e-4*xkahan[i-5]*xkahan[i-1]**4,\
                  + 0.2484e-4*xkahan[i-3]**2*xkahan[i-2]**3,\
                  + 0.1238e-2*xkahan[i-2]*xkahan[i-1]**2,\
                  + 0.4353e-4*xkahan[i-5]**4,\
                  + 0.2258e-2*xkahan[i-5]*xkahan[i-2]*xkahan[i-1]**2,\
                  + 0.3123e-4*xkahan[i-4]**5,\
                  + 0.7531e-3*xkahan[i-1]**4,\
                  - 0.2703e-2*xkahan[i-3]**2*xkahan[i-1]**2,\
                  - 0.7807e-3*xkahan[i-1]**3,\
                  - 0.7077e-4*xkahan[i-3]**2*xkahan[i-2]**2*xkahan[i-1],\
                  - 0.3304e-3*xkahan[i-3]*xkahan[i-2]**3,\
                  - 0.8847e-2*xkahan[i-5]*xkahan[i-1],\
                  + 0.7631e-2*xkahan[i-4]*xkahan[i-1],\
                  - 0.3870e-4*xkahan[i-5]**3*xkahan[i-1]**2,\
                  + 0.4676e-3*xkahan[i-3]**3*xkahan[i-1]]) )



#==============================================================================
#===================== EXPOENTE DE LYAPUNOV ===================================
#---------------------------- + NUMPY +----------------------------------------


# PARÂMETROS PARA PONTOS DA LINEARREGRESSION E LBE
start=5 #ponto que começa calcular o LBE (de preferência o mesmo de h)
vec1=0 #primeiro ponto da reta
vec2=1800#último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEpy=abs(xref[start:N]-xnumpy[start:N])/2 #calcula o LBE
LBEpy=np.array(LBEpy,dtype=np.float64) #transformando para array tipo np
LBEnumpy=np.log2(LBEpy) #gera o LBE no gráfico de log2

LLEnumpy=LinearRegression()
LLEnumpy.fit(vecN.reshape(-1,1),LBEnumpy[vec1:vec2])
LLEpy_coef=LLEnumpy.coef_
LLEpy_inter=LLEnumpy.intercept_

Ncnumpy=int(abs(np.log2(18.304/2)+np.log2(10**16))/LLEpy_coef)

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
plt.text(Ncnumpy+1, -59, 'N$_c$= %s'%Ncnumpy,fontsize=16, color="red" )
plt.text(500, min(LBEnumpy), '%0.4fn'%LLEpy_coef,fontsize=17 )
plt.text(1155, min(LBEnumpy), '%0.2f'%LLEpy_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEnumpy)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{numpy}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('rossler_lbe_numpy.pdf', format='pdf')

#---------------------------- + NUMPY KAHAN +-----------------------------------

vec1=0 #primeiro ponto da reta
vec2=2300 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEkhpy=abs(xref[start:N]-xkahan[start:N])/2 #calcula o LBE
LBEkhpy=np.array(LBEkhpy,dtype=np.float64) #transformando para array tipo np
LBEkahanpy=np.log2(LBEkhpy) #gera o LBE no gráfico de log2

LLEkahanpy=LinearRegression()
LLEkahanpy.fit(vecN.reshape(-1,1),LBEkahanpy[vec1:vec2])
LLEkhpy_coef=LLEkahanpy.coef_
LLEkhpy_inter=LLEkahanpy.intercept_


Nckahanpy=int(abs(np.log2(18.304/2)+np.log2(10**16))/LLEkahanpy.coef_)

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
plt.text(Nckahanpy+1, -59.5, 'N$_c$= %s'%Nckahanpy,fontsize=16, color="red" )
plt.text(500, min(LBEkahanpy)+0.5, '%0.4fn'%LLEkhpy_coef,fontsize=17 )
plt.text(1155, min(LBEkahanpy)+0.5, '%0.2f'%LLEkhpy_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEkahanpy)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{kahan}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('rossler_lbe_pykahan.pdf', format='pdf')

plt.figure(figsize=(8, 6))
plt.plot(xref,color='black',linewidth=2, alpha=0.8,marker = None)
plt.plot(xnumpy,color='red',linewidth=1,marker=None)
plt.axvline(x=Ncnumpy, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Ncnumpy+3,-11, 'N$_c$= %s'%Ncnumpy,fontsize=16 )
plt.axvline(x=2050, color='black',linestyle="solid")
plt.axhline(y=-12, color='black',linestyle="solid")
plt.axis([2050, 2250, -12, 15])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{numpy}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('rossler_ref_numpy.pdf', format='pdf')

plt.figure(figsize=(8, 6))
plt.plot(xref,color='black',linewidth=2, alpha=0.8,marker = None)
plt.plot(xkahan,color='red',linewidth=1,marker=None)
plt.axvline(x=Nckahanpy, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Nckahanpy+3,-11, 'N$_c$= %s'%Nckahanpy,fontsize=16 )
plt.axvline(x=2300, color='black',linestyle="solid")
plt.axhline(y=-12, color='black',linestyle="solid")
plt.axis([2300, 2600, -12, 15])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{kahan}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('rossler_ref_pykahan.pdf', format='pdf')




#==============================================================================
#=======================-----NRMSE-------======================================
#==============================================================================
n=int(Ncnumpy/10)
#----------------------------python--------------------------------------------
NRMSE_numpy=[]
for i in range(1,11):
    NRMSE_numpy=np.append(NRMSE_numpy, np.sqrt(np.sum((xref[0:i*n] - xnumpy[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xnumpy[0:i*n])**2)))))


NRMSE_kahanpy=[]
for i in range(1,11):
    NRMSE_kahanpy=np.append(NRMSE_kahanpy, np.sqrt(np.sum((xref[0:i*n] - xkahan[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xkahan[0:i*n])**2)))))


print('============================================================')
print('PYTHON POR PARTES-----> NRMSE REFERÊNCIA X NUMPY')
NRMSE_numpy=np.array(NRMSE_numpy.tolist(),dtype=np.float64) #transformar mpf to numpy array
NRMSE_kahanpy=np.array(NRMSE_kahanpy.tolist(),dtype=np.float64) #transformar mpf to numpy array
title=[" ","NUMPY","PYTHON KAHAN"]
it=np.linspace(n,Nckahanpy, num=10)
it=it.astype(int)
#NRMSE=[[NRMSE_numpy] [NRMSE_kahanpy]]
NRMSE_2=pd.DataFrame({ 'numpy': NRMSE_numpy,
                    'python kahan': NRMSE_kahanpy})
NRMSE_2.style
NRMSE_2=NRMSE_2.set_index([pd.Index(it)])
#NRMSE.style
print(NRMSE_2)





#==============================================================================
#===================== EXPOENTE DE LYAPUNOV ===================================
#---------------------------- + MATALB +---------------------------------------


# PARÂMETROS PARA PONTOS DA LINEARREGRESSION E LBE
start=5 #ponto que começa calcular o LBE (de preferência o mesmo de h)
vec1=0 #primeiro ponto da reta
vec2=2000 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEmat=abs(xref[start:N]-xrossler_normal[start:N])/2 #calcula o LBE
LBEmat=np.array(LBEmat,dtype=np.float64) #transformando para array tipo np
LBEmatlab=np.log2(LBEmat) #gera o LBE no gráfico de log2

LLEmatlab=LinearRegression()
LLEmatlab.fit(vecN.reshape(-1,1),LBEmatlab[vec1:vec2])
LLEmat_coef=LLEmatlab.coef_
LLEmat_inter=LLEmatlab.intercept_

Ncmatlab=int(abs(np.log2(1.284/2)+np.log2(10**16))/LLEmat_coef)

print("-----------------------")
print("LLE Matlab")
print(LLEmat_coef)
print("Tempo Matlab")
print(Ncmatlab)

# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEmatlab,color='black',linewidth=1)
plt.plot(vecN,(vecN*LLEmat_coef+LLEpy_inter),color='red',linewidth=3)
plt.axvline(x=Ncmatlab, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEmatlab)-3, color='black',linestyle="solid")
plt.text(Ncmatlab-4, -59, 'N$_c$= %s'%Ncmatlab,fontsize=16, color="red" )
plt.text(500, min(LBEmatlab)+1, '%0.4fn'%LLEmat_coef,fontsize=17 )
plt.text(1200, min(LBEmatlab)+1, '%0.2f'%LLEmat_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEmatlab)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{matlab}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('rossler_lbe_matlab.pdf', format='pdf')
#---------------------------- + MATLAB KAHAN +-----------------------------------

vec1=0 #primeiro ponto da reta
vec2=2200 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression


LBEkhmat=abs(xref[start:N]-xrossler_kahan[start:N])/2 #calcula o LBE
LBEkhmat=np.array(LBEkhmat,dtype=np.float64) #transformando para array tipo np
LBEkahanmat=np.log2(LBEkhmat) #gera o LBE no gráfico de log2

LLEkahanmat=LinearRegression()
LLEkahanmat.fit(vecN.reshape(-1,1),LBEkahanmat[vec1:vec2])
LLEkhmat_coef=LLEkahanmat.coef_
LLEkhmat_inter=LLEkahanmat.intercept_


Nckahanmat=int(abs(np.log2(1.119/2)+np.log2(10**16))/LLEkahanmat.coef_)

print("-----------------------")
print("LLE matlab KAHAN")
print(LLEkhmat_coef)
print("Tempo matlab KAHAN")
print(Nckahanmat)



# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEkahanmat,color='black',linewidth=1)
plt.plot(vecN,(vecN*LLEkhmat_coef+LLEkhmat_inter),color='red',linewidth=2)
plt.axvline(x=Nckahanmat, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEkahanmat)-3, color='black',linestyle="solid")
plt.text(Nckahanmat-4, -59, 'N$_c$= %s'%Nckahanmat,fontsize=16, color="red" )
plt.text(500, min(LBEkahanmat)+1, '%0.4fn'%LLEkhmat_coef,fontsize=17 )
plt.text(1200, min(LBEkahanmat)+1, '%0.2f'%LLEkhmat_inter,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEkahanmat)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{kahan}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('rossler_lbe_matkahan.pdf', format='pdf')

#==============================================================================
#=======================-----NRMSE-------======================================
#==============================================================================
n=int(Nckahan/10)
#----------------------------matlab--------------------------------------------
#matlab
NRMSE_normal=[]
for i in range(1,11):
    NRMSE_normal=np.append(NRMSE_normal, np.sqrt(np.sum((xref[0:i*n] - xrossler_normal[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xrossler_normal[0:i*n])**2)))))


NRMSE_kahan=[]
for i in range(1,11):
    NRMSE_kahan=np.append(NRMSE_kahan, np.sqrt(np.sum((xref[0:i*n] - xrossler_kahan[0:i*n])**2))/(np.sqrt(np.sum((xref[0:i*n] - np.mean(xrossler_kahan[0:i*n])**2)))))



print('============================================================')
print('MATALAB POR PARTES-----> NRMSE REFERÊNCIA X COM KAHAN')
NRMSE_kahan=np.array(NRMSE_kahan.tolist(),dtype=np.float64) #transformar mpf to numpy array
NRMSE_normal=np.array(NRMSE_normal.tolist(),dtype=np.float64) #transformar mpf to numpy array
title=[" ","MATLAB","MATLAB KAHAN"]
it=np.linspace(n,Nckahanmat, num=10)
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
plt.plot(xrossler_normal,color='red',linewidth=1,marker=None)
plt.axvline(x=Ncmatlab, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Ncmatlab+3,-11, 'N$_c$= %s'%Ncmatlab,fontsize=16 )
plt.axvline(x=2050, color='black',linestyle="solid")
plt.axhline(y=-12, color='black',linestyle="solid")
plt.axis([2050, 2300, -12, 15])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N ',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{matlab}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('rossler_ref_matlab.pdf', format='pdf')

plt.figure(figsize=(8, 6))
plt.plot(xref,color='black',linewidth=2, alpha=0.8,marker = None)
plt.plot(xrossler_kahan,color='red',linewidth=1,marker=None)
plt.axvline(x=Nckahanmat, color='green',linestyle='--',alpha=0.7,linewidth=2)
plt.text(Nckahanmat+3,-11, 'N$_c$= %s'%Nckahanmat,fontsize=16 )
plt.axvline(x=2150, color='black',linestyle="solid")
plt.axhline(y=-12, color='black',linestyle="solid")
plt.axis([2150, 2500, -12, 15])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('x$_{ref}$ , x$_{kahan}$',fontsize=17)
plt.grid(True, linestyle=':')
plt.box(False)
plt.savefig('rossler_ref_matkahan.pdf', format='pdf')


