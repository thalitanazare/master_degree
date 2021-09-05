%======Routine to calculate NRMSE Jerk Circuit in LTspice 4 and  17 =======
%Author:Thalita Emanuelle de Nazaré
%Year:2020 
%Reproducibility of the chaotic circuit
%------------------------Computadores Utilizados---------------------------
%Computer 4
%--------------------------------------------------------------------------
clc;clear all;close all;
format long
%-----------------Loading the vector exported from LTspice-----------------
%The txt files must be in the same folder as the code

load Jerk4.txt 
load Jerk17.txt 

%Voltage Vectors
z4 = Jerk4(:,2);
z17 = Jerk17(:,2); 

%Time Vectors
t4 = Jerk4(:,1); 
t17 = Jerk17(:,1); 

%-----------------------------Interpolation--------------------------------
temp=linspace(0,0.1,22900);
z4temp=interp1(t4,z4,temp)';
z17temp=interp1(t17,z17,temp)';
%==========================================================================
%-----------------------      NMRSE - Computer     ------------------------
%Jerk4 and Jerk17
numc1=(z17temp-z4temp)'*(z17temp-z4temp); %calculating the numerator
denc1=(z17temp-mean(z4temp))'*(z17temp-mean(z4temp)); %calc the denominator
nrmsec1=sqrt(numc1)/sqrt(denc1); %calc the nrmse
disp('NRMSE LTspice 4 e 17')
disp(nrmsec1')
%------------------for different simulation intervals----------------------
%Computer 1 
for k=1:10
    NUM1=(z17temp(1:k*2290)-z4temp(1:k*2290))'*(z17temp(1:k*2290)-z4temp(1:k*2290)) ;
    DEN1=(z17temp(1:k*2290)-mean(z4temp(1:k*2290)))'*(z17temp(1:k*2290)-mean(z4temp(1:k*2290))); 
    NRMSE1(k)=sqrt(NUM1)/sqrt(DEN1); 
end
disp('NRMSE LTspice 4 e 17 - Partes')
disp(NRMSE1')



%==========================================================================


%==========================================================================
%FIGURES
figure(1)
subplot(3,1,1)
plot(t4,z4,'-','LineWidth',2,'Color',[1 0 0]);
hold on
plot(t17,z17,'-','LineWidth',1,'Color',[0 0 0]);
xlabel('Tempo (s)','FontSize',20,'FontName','Times');
ylabel('Tensão (V)','FontSize',20,'FontName','Times');
set(gca,'fontsize',20,'FontName','Times')
grid on;
box off;

subplot(3,1,2)
plot(t4,z4,'.','LineWidth',1,'Color',[1 0 0])
hold on
plot(t17,z17,'.','LineWidth',1,'Color',[0 0 0])
xlabel('Tempo (s)','FontSize',20,'FontName','Times');
ylabel('Tensão (V)','FontSize',20,'FontName','Times');
set(gca,'fontsize',20,'FontName','Times')
grid on;
box off;

subplot(3,1,3)
plot(t4,z4,'.','LineWidth',1,'Color',[1 0 0])
hold on
plot(t17,z17,'.','LineWidth',1,'Color',[0 0 0])
xlabel('Tempo (s)','FontSize',20,'FontName','Times');
ylabel('Tensão (V)','FontSize',20,'FontName','Times');
set(gca,'fontsize',20,'FontName','Times')
grid on;
box off;
%--------------------------------------------------------------------------