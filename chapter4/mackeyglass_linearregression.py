# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 15:33:30 2021

@author: thali
"""

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
import pandas as pd
import matplotlib.pyplot as plt
import math as mt
from decimal import *
plt.rcParams["font.family"] = "Times New Roman" 
import mpmath as mp
from mpmath import *
from multivar_horner.multivar_horner import HornerMultivarPolynomial
from sympy import Symbol, solveset, S
import random as r
import scipy.stats as stats
from sklearn.linear_model import LinearRegression

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

# PRECISÃO: PARA SISTEMA DE 128 bits
mp.dps=34 # casas decimais -> 116 bits

# VETORES DE ENTRADA
cond=[] #condição inicial
LBEnum=[] #para acumulo dos vetores de LBEnumpy
Tsnum=[] #para acumulo dos valores de tempo crítico
LLEnum=[] #para acumulo dos valores de inclinação da reta: LLE
LLEnumintercept=[] #para acumulo dos valores de intercept
LBEka=[]
Tska=[]
LLEka=[]
LLEkaintercept=[]
xref1=[]
xnumpy1=[]
xkahan1=[]

# PARÂMETROS PARA PONTOS DA LINEARREGRESSION E LBE
start=10 #ponto que começa calcular o LBE (de preferência o mesmo de h)
vec1=30 #primeiro ponto da reta
vec2=2200 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression

vec3=0 #primeiro ponto da reta
vec4=2000 #último ponto da reta
num2=vec4-vec3 #tamanho do vetor
veck=np.linspace(vec3,vec4,num=num2) #vetor de iterações para a linearregression

# PARÂMETROS PARA SIMULAÇÃO
h=10 #número de condições iniciais necessárias para começar a simular (número de atrasos)
N=5000 #número de iterações
r.seed(3) #valor que determina o ponto de números aleatórios: permite sempre sortear o mesmo conjunto

# INÍCIO SIMULAÇÃO: Gera condições Iniciais
for i in range(1,301): #pq preciso de 10 cond inciais e rodar 30 vezes
    cond=np.append(cond, r.uniform(0.25,0.35))#0.3 0.301 funcionou para 5 vezes
    
# Simula por 'x' vezes, ou seja, faz 'x' testes
for j in range(0,15):
    x0=[cond[10*j], cond[10*j+1], cond[10*j+2], cond[10*j+3], cond[10*j+4],cond[10*j+5], cond[10*j+6], cond[10*j+7], cond[10*j+8], cond[10*j+9]]
    xref=[mp.mpf(x0[0]),mp.mpf(x0[1]),mp.mpf(x0[2]),mp.mpf(x0[3]),mp.mpf(x0[4]),mp.mpf(x0[5]),mp.mpf(x0[6]),mp.mpf(x0[7]),mp.mpf(x0[8]),mp.mpf(x0[9])]
    xnumpy=x0
    xkahan=x0
#==============================================================================
#=======================....SIMULAÇÃO....======================================
#==============================================================================
        
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

    xref1.append(xref)
    xnumpy1.append(xnumpy)
    xkahan1.append(xkahan)
#==============================================================================
#===================== EXPOENTE DE LYAPUNOV ===================================
#---------------------------- + NUMPY +----------------------------------------

    # calculando lower bound error
    LBEpy=abs(xref[start:N]-xnumpy[start:N])/2 #calcula o LBE
    LBEpy=np.array(LBEpy,dtype=np.float) #transformando para array tipo np
    LBEnumpy=np.log2(LBEpy) #gera o LBE no gráfico de log2
    LBEnum.append(LBEnumpy) #acumula para média do gráfico
    # calculando o expoente de lyapunov
    LLEnumpy=LinearRegression()
    LLEnumpy.fit(vecN.reshape(-1,1),LBEnumpy[vec1:vec2]) #calcula a reta
    LLEnum=np.append(LLEnum,LLEnumpy.coef_) #acumula para média
    LLEnumintercept=np.append(LLEnumintercept, LLEnumpy.intercept_)
    # calculando o tempo crítico
    Tsnumpy=abs(np.log10(1.486/2)+np.log2(10**16))/LLEnumpy.coef_
    Tsnum=np.append(Tsnum,Tsnumpy)

#==============================================================================
#===================== EXPOENTE DE LYAPUNOV ===================================
#---------------------------- + KAHAN +----------------------------------------    

    # calculando lower bound error
    LBEk=abs(xref[start:N]-xkahan[start:N])/2 #calcula o LBE
    LBEk=np.array(LBEk,dtype=np.float) #transformando para array tipo np
    LBEkahan=np.log2(LBEk) #gera o LBE no gráfico de log2
    LBEka.append(LBEkahan) #acumula para média do gráfico
    # calculando o expoente de lyapunov
    LLEkahan=LinearRegression()
    LLEkahan.fit(veck.reshape(-1,1),LBEkahan[vec3:vec4]) #calcula a reta
    LLEka=np.append(LLEka,LLEkahan.coef_) #acumula para média
    LLEkaintercept=np.append(LLEkaintercept, LLEkahan.intercept_)
    # calculando o tempo crítico
    Tskahan=abs(np.log10(1.486/2)+np.log2(10**16))/LLEkahan.coef_
    Tska=np.append(Tska,Tskahan)


#==============================================================================
#======================== ++ VALORES FINAIS ++ ================================


# Média de todos os valores para mostrar o resultado final
LBEnumpy=np.mean(LBEnum,axis=0) 
LLEnumpy=np.mean(LLEnum)
LLEnumintercept=np.mean(LLEnumintercept)
Tsnumpy=int(np.mean(Tsnum))
print("-------------MÉDIA----------")
print("LLE numpy")
print(LLEnumpy)
print("Tempo numpy")
print(int(Tsnumpy))

LBEkahan=np.mean(LBEka,axis=0) 
LLEkahan=np.mean(LLEka)
LLEkaintercept=np.mean(LLEkaintercept)
Tskahan=int(np.mean(Tska))
print("-------------MÉDIA----------")
print("LLE numpy")
print(LLEkahan)
print("Tempo numpy")
print(int(Tskahan))



# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEnumpy,color='black',linewidth=0.5)
plt.plot(vecN,(vecN*LLEnumpy+LLEnumintercept),color='red',linewidth=3)
plt.axvline(x=Tsnumpy, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEnumpy)-3, color='black',linestyle="solid")
plt.text(Tsnumpy+1, -57, 'N$_c$= %s'%Tsnumpy,fontsize=16, color="red" )
plt.text(570, min(LBEnumpy)-2, '%0.3fn'%LLEnumpy,fontsize=17 )
plt.text(1150, min(LBEnumpy)-2, '%0.3f'%LLEnumintercept,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEnumpy)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{numpy}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('mg_media_numpy.pdf', format='pdf')

# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEkahan,color='black',linewidth=0.5)
plt.plot(veck,(veck*LLEkahan+LLEkaintercept),color='red',linewidth=3)
plt.axvline(x=Tskahan, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEkahan)-3, color='black',linestyle="solid")
plt.text(Tskahan+1, -58, 'N$_c$= %s'%Tskahan,fontsize=16, color='red' )
plt.text(570, min(LBEkahan)-2, '%0.3fn'%LLEkahan,fontsize=17 )
plt.text(1150, min(LBEkahan)-2, '%0.3f'%LLEkaintercept,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEkahan)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{kahan}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('mg_media_kahan.pdf', format='pdf')
























