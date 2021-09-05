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
format shortg
intvalinit('displayinfsup')

p1=0.003;
p2=0.003;
p3=0.003;


A=[infsup(0.5,0.5+p1) infsup(0,0) intval(0.5);
   intval(1) infsup(0.5,0.5+p2) intval(1);
   intval(0.5) intval(0) infsup(0.25,0.25+p3)];

C=intval(eye(size(A,2)));

I=eye(size(A,1)^2);

G=I-kron(A',A');

c=vec(C);
x=verifylss(G,c);

X3=[x(1) x(2) x(3);
    x(4) x(5) x(6);
    x(7) x(8) x(9)];

X2=[x(1) x(2);
    x(4) x(5)];

X1=x(1);

%----------------------------------------
disp('Matriz X - Considerando Incertezas Numéricas')
disp(X3)

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
n=size(A,2);
x=intval(randn(n,1));
t=100;
for k=1:t
    x(:,k+1)=A*x(:,k);
    Lyap(:,k)=x(:,k)'*X3*x(:,k);

end
figure()
for i=1:n
 T=1:size(x,2);
plot(T,x(i,:));
hold on
end

figure()
plot(Lyap)
% for i=1:1:10 
% n=size(A,2);
% x=intval(randn(n,1));
% t=100;
% for k=1:t
%     x(:,k+1)=A*x(:,k);
%     V(:,k)=x(:,k)'*X3*x(:,k);
% end
% k=1:1:101; 
% plot(k,x)
% figure
% plot(V)
% pause
% close
% clear x
% end