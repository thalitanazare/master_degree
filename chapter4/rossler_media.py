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

# PARÂMETROS PARA PONTOS DA LINEARREGRESSION E LBE
start=3 #ponto que começa calcular o LBE (de preferência o mesmo de h)
vec1=0 #primeiro ponto da reta
vec2=2200 #último ponto da reta
num1=vec2-vec1 #tamanho do vetor
vecN=np.linspace(vec1,vec2,num=num1) #vetor de iterações para a linearregression

vec3=100 #primeiro ponto da reta
vec4=2400 #último ponto da reta
num2=vec4-vec3 #tamanho do vetor
veck=np.linspace(vec3,vec4,num=num2) #vetor de iterações para a linearregression

# PARÂMETROS PARA SIMULAÇÃO
h=5 #número de condições iniciais necessárias para começar a simular (número de atrasos)
N=5000 #número de iterações
r.seed(3) #valor que determina o ponto de números aleatórios: permite sempre sortear o mesmo conjunto


for i in range(1,151): #preciso de 5 cond iniciais e rodar 30 vezes
    cond=np.append(cond, r.uniform(0.9,1.1))

for j in range(0,10):
    x0=[cond[5*j], cond[5*j+1], cond[5*j+2], cond[5*j+3], cond[5*j+4]]
    xref=[mp.mpf(x0[0]),mp.mpf(x0[1]),mp.mpf(x0[2]),mp.mpf(x0[3]),mp.mpf(x0[4])]

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



      #linhas 75 77 79 83

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

#EXPOENTE 
    start=5

#    NUMPY
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
    Tsnumpy=abs(np.log10(6.22/2)+np.log2(10**16))/LLEnumpy.coef_
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
    Tskahan=abs(np.log10(4.66/2)+np.log2(10**16))/LLEkahan.coef_
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
plt.plot(LBEnumpy,color='black',linewidth=1)
plt.plot(vecN,(vecN*LLEnumpy+LLEnumintercept),color='red',linewidth=3)
plt.axvline(x=Tsnumpy, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEnumpy)-3, color='black',linestyle="solid")
plt.text(Tsnumpy+1, -55, 'N$_c$= %s'%Tsnumpy,fontsize=16, color="red" )
plt.text(500, min(LBEnumpy)+5, '%0.4fn'%LLEnumpy,fontsize=17 )
plt.text(1150, min(LBEnumpy)+5, '%0.2f'%LLEnumintercept,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEnumpy)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{numpy}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('rossler_media_numpy.pdf', format='pdf')

# Figura: considerando já os valores com média
plt.figure(figsize=(8, 6))
plt.plot(LBEkahan,color='black',linewidth=1)
plt.plot(veck,(veck*LLEkahan+LLEkaintercept),color='red',linewidth=3)
plt.axvline(x=Tskahan, color='grey',linestyle='--',alpha=0.6,linewidth=2)
plt.axvline(x=0, color='black',linestyle="solid")
plt.axhline(y=min(LBEkahan)-3, color='black',linestyle="solid")
plt.text(Tskahan+1, -55, 'N$_c$= %s'%Tskahan,fontsize=16, color='red' )
plt.text(500, min(LBEkahan)+5, '%0.4fn'%LLEkahan,fontsize=17 )
plt.text(1150, min(LBEkahan)+5, '%0.2f'%LLEkaintercept,fontsize=17 )
#plt.text(1100, min(graf_LBEN)+5, '%0.3f'%media_LLEparte ,fontsize=17 )
#plt.vlines(0, 4000, -20, linestyles ="solid", colors ="k")
plt.axis([0, N, min(LBEkahan)-3, 2])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('N',fontsize=17)
plt.ylabel('log$_2$(|x$_{ref}$ - x$_{kahan}$|)',fontsize=17)
plt.grid(True,alpha=0.8,linestyle="--")
plt.box(False)
plt.savefig('rossler_media_kahan.pdf', format='pdf')