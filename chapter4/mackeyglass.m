% Mackey-Glass com e sem kahan
% Lembrar de substituir os dados lá na pasta análise
clear all
clc
close all
format long
x=[0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01];
xk=[0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01];
s=[0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01];
N=6000;
for k=11:N
x(k)=0.24662*10*x(k-1) - 0.16423*10*x(k-2) + 0.60992*x(k-3) + 0.73012e-1*x(k-5)^2*x(k-10)^2 + 0.38566*x(k-3)*x(k-10) + 0.66999*x(k-1)*x(k-10)^2+ 0.88364*x(k-1)^3- 0.67300*x(k-4)*x(k-10)^2- 0.11929*10*x(k-1)^2 - 0.50451e-1*x(k-4)*x(k-5) - 0.24765*x(k-1)^4 + 0.42081*x(k-4)*x(k-9)*x(k-10)^2- 0.70406*x(k-1)*x(k-10)^3- 0.14089*x(k-5)*x(k-8)^2+ 0.14807*x(k-1)*x(k-7)*x(k-10);
end

for i=11:N
xk=[0.24662*10*s(i-1); - 0.16423*10*s(i-2) ; + 0.60992*s(i-3) ; + 0.73012e-1*s(i-5)^2*s(i-10)^2 ; + 0.38566*s(i-3)*s(i-10) ; + 0.66999*s(i-1)*s(i-10)^2; + 0.88364*s(i-1)^3;- 0.67300*s(i-4)*s(i-10)^2;- 0.11929*10*s(i-1)^2 ;- 0.50451e-1*s(i-4)*s(i-5); - 0.24765*s(i-1)^4 ; + 0.42081*s(i-4)*s(i-9)*s(i-10)^2;- 0.70406*s(i-1)*s(i-10)^3;- 0.14089*s(i-5)*s(i-8)^2; + 0.14807*s(i-1)*s(i-7)*s(i-10)];
s(i)=ksum(xk);
end

plot(x,'k')
hold on
plot(s)

filex = fopen('xmackeyglass_normal.txt','w');
fprintf(filex,'%12.15f\n',x);
fclose(filex);

filexk = fopen('xmackeyglass_kahan.txt','w');
fprintf(filexk,'%12.15f\n',s);
fclose(filexk);