% ==================================================================== %
%                       COMPUTAÇÃO ARITMÉTICA
%                 Análise da Estabilidade Intervalar                   %
% -------------------------------------------------------------------- %
% MÉTODO PAULO THALITA
%EXEMPLO SBAI 2015 HELOISE
% -------------------------------------------------------------------- %
% Alunos: Thalita Emanuelle e Paulo 
% Data: 01/12/2019
% ==================================================================== %
clear all
close all
clc
format long
intvalinit('displayinfsup')
m1=2.77; %considera como constante
m2=2.59;

p=0.006;
c11=1.2;
c22=0.2;
k11=200;
k22=390;
k33=30;

T=0.0001;
m1=infsup(m1-m1*p,m1+m1*p);
m2=infsup(m2-m2*p,m2+m2*p);
c1=infsup(c11-c11*p,c11+c11*p);
c2=infsup(c22-c22*p,c22+c22*p);
k1=infsup(k11-k11*p,k11+k11*p);
k2=infsup(k22-k22*p,k22+k22*p);
k3=infsup(k33-k33*p,k33+k33*p);

A=[intval(0) intval(1) intval(0) intval(0);
    -T*(k1+k2)/m1 1-T*c1/m1 T*k2/m1 intval(0);
    intval(0) intval(0) intval(0) intval(1);
    T*k2/m2 intval(0) -T*(k2+k3)/m2 1-T*c2/m2];
C=intval(eye(size(A,2)));

I=eye(size(A,1)^2);

G=I-kron(A',A');


c=vec(C);
x=verifylss(G,c);
x11=x(1);
x12=x(2);
x13=x(3);
x14=x(4);

x21=x(5);
x22=x(6);
x23=x(7);
x24=x(8);

x31=x(9);
x32=x(10);
x33=x(11);
x34=x(12);

x41=x(13);
x42=x(14);
x43=x(15);
x44=x(16);



X4=[x11 x12 x13 x14;
    x21 x22 x23 x24;
    x31 x32 x33 x34;
    x41 x42 x43 x44];
X3=[x11 x12 x13;
    x21 x22 x23;
    x31 x32 x33];
X2=[x11 x12;
    x21 x22];
X1=x11;

%----------------------------------------
disp('Matriz X - Considerando Incertezas Numéricas')
disp(X4)

disp('----------------------------------')
disp('-----------Determinante 4--------')
det4=calculo_determinante.calc_det(X4);
disp(det4)

disp('----------------------------------')
disp('-----------Determinante 3--------')
det3=calculo_determinante.calc_det(X3);
disp(det3)

disp('----------------------------------')
disp('-----------Determinante 2--------')
det2=calculo_determinante.calc_det(X2);
disp(det2)

disp('----------------------------------')
disp('-----------Determinante 1--------')
disp(X1)

%% ===============================================================

