% Mapa de henon com e sem kahan

clear all
clc
close all
format long
x=[0.01;0.01;0.01];
xk=[0.01;0.01;0.01];
s=[0.01;0.01;0.01];
N=300;
for k=4:N
x(k)=-1.3772*x(k-1)^2 + 0.96958 + 0.00083*x(k-1)^2*x(k-3) + 0.3410*x(k-2) - 0.03352*x(k-1)*x(k-2)*x(k-3)- 0.04836*x(k-1)*x(k-3)^2;
end

for i=4:N
xk=[-1.3772*s(i-1)^2;0.96958;0.00083*s(i-1)^2*s(i-3);0.3410*s(i-2);- 0.03352*s(i-1)*s(i-2)*s(i-3);- 0.04836*s(i-1)*s(i-3)^2];
s(i)=ksum(xk);
%sum(k)=ksum(xk);
end

plot(x,'k')
hold on
plot(s)

filex = fopen('xhenon_normal.txt','w');
fprintf(filex,'%12.15f\n',x);
fclose(filex);

filexk = fopen('xhenon_kahan.txt','w');
fprintf(filexk,'%12.15f\n',s);
fclose(filexk);
