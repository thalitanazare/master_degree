%========Routine to calculate NRMSE Jerk Circuit with LTspiceXVII==========
%Author:Thalita Emanuelle de Nazare
%Year:2020 
%Reproducibility of the chaotic circuit 
%------------------------Computadores Utilizados---------------------------
%Computer 1 - main core i5
%Computer 2 - Lecom dual core 2.5Ghz
%Computer 3 - Lecom core i5-3570
%Computer 4 - i5-4210U
%--------------------------------------------------------------------------
clc;clear all;close all;
format long
%-----------------Loading the vector exported from LTspice------------------
%The txt files must be in the same folder as the code

load Jerk_1cond.txt 
load Jerk_2cond.txt 
load Jerk_3cond.txt 
load Jerk_4cond.txt 
load Jerk_circ.txt %voltage data collected in the real circuit.
%load DIODO.txt
%Voltage Vectors
z1 = Jerk_1cond(:,2);
z2 = Jerk_2cond(:,2); 
z3 = Jerk_3cond(:,2);
z4 = Jerk_4cond(:,2);
zC = Jerk_circ(:,1);
%Time Vectors
t1 = Jerk_1cond(:,1); 
t2 = Jerk_2cond(:,1); 
t3 = Jerk_3cond(:,1);
t4 = Jerk_4cond(:,1);
tC = linspace(0,0.1,10000);
%-----------------------------Interpolation--------------------------------
temp=linspace(0,0.1,4700);
z1temp=interp1(t1,z1,temp)';
z2temp=interp1(t2,z2,temp)';
z3temp=interp1(t3,z3,temp)';
z4temp=interp1(t4,z4,temp)';
zCtemp=interp1(tC,zC,temp)';
%==========================================================================
%-----------------------      NMRSE - Computer     ------------------------
%Computer 1
numc1=(zCtemp-z1temp)'*(zCtemp-z1temp); %calculating the numerator
denc1=(zCtemp-mean(z1temp))'*(zCtemp-mean(z1temp)); %calc the denominator
nrmsec1=sqrt(numc1)/sqrt(denc1); %calc the nrmse
%Computer 2
numc2=(zCtemp-z2temp)'*(zCtemp-z2temp); 
denc2=(zCtemp-mean(z2temp))'*(zCtemp-mean(z2temp)); 
nrmsec2=sqrt(numc2)/sqrt(denc2) ;
%Computer 3
numc3=(zCtemp-z3temp)'*(zCtemp-z3temp); 
denc3=(zCtemp-mean(z3temp))'*(zCtemp-mean(z3temp)); 
nrmsec3=sqrt(numc3)/sqrt(denc3) ;
%Computer 4
numc4=(zCtemp-z4temp)'*(zCtemp-z4temp); 
denc4=(zCtemp-mean(z4temp))'*(zCtemp-mean(z4temp)); 
nrmsec4=sqrt(numc4)/sqrt(denc4) ;
%------------------for different simulation intervals----------------------
%Computer 1 
for k=1:10
    NUM1=(zCtemp(1:k*470)-z1temp(1:k*470))'*(zCtemp(1:k*470)-z1temp(1:k*470)) ;
    DEN1=(zCtemp(1:k*470)-mean(z1temp(1:k*470)))'*(zCtemp(1:k*470)-mean(z1temp(1:k*470))); 
    NRMSE1(k)=sqrt(NUM1)/sqrt(DEN1); 
end
%Computer 2 
for k=1:10
    NUM2=(zCtemp(1:k*470)-z2temp(1:k*470))'*(zCtemp(1:k*470)-z2temp(1:k*470));
    DEN2=(zCtemp(1:k*470)-mean(z2temp(1:k*470)))'*(zCtemp(1:k*470)-mean(z2temp(1:k*470)));
    NRMSE2(k)=sqrt(NUM2)/sqrt(DEN2); 
end
%Computer 3 
for k=1:10
    NUM3=(zCtemp(1:k*470)-z3temp(1:k*470))'*(zCtemp(1:k*470)-z3temp(1:k*470));
    DEN3=(zCtemp(1:k*470)-mean(z3temp(1:k*470)))'*(zCtemp(1:k*470)-mean(z3temp(1:k*470)));
    NRMSE3(k)=sqrt(NUM3)/sqrt(DEN3); 
end

%Computer 4 
for k=1:10
    NUM4=(zCtemp(1:k*470)-z4temp(1:k*470))'*(zCtemp(1:k*470)-z4temp(1:k*470));
    DEN4=(zCtemp(1:k*470)-mean(z4temp(1:k*470)))'*(zCtemp(1:k*470)-mean(z4temp(1:k*470)));
    NRMSE4(k)=sqrt(NUM4)/sqrt(DEN4);
end


%==========================================================================
%------------------  To Calculate NMRSE Between Computers  ----------------
%Computer 4 with 1
num1=(z4temp-z1temp)'*(z4temp-z1temp); 
den1=(z4temp-mean(z1temp))'*(z4temp-mean(z1temp)); 
nrmse1=sqrt(num1)/sqrt(den1);

%Computer 4 with 2
num2=(z4temp-z2temp)'*(z4temp-z2temp) ;
den2=(z4temp-mean(z2temp))'*(z4temp-mean(z2temp)); 
nrmse2=sqrt(num2)/sqrt(den2);

%Computer 4 with 3
num3=(z4temp-z3temp)'*(z4temp-z3temp); 
den3=(z4temp-mean(z3temp))'*(z4temp-mean(z3temp)); 
nrmse3=sqrt(num3)/sqrt(den3);

%==========================================================================
%------------------  To Calculate NMRSE Between Computers  ----------------
%Computer 4 with 1
for k=1:10
num1=(z4temp(1:k*470)-z1temp(1:k*470))'*(z4temp(1:k*470)-z1temp(1:k*470)); 
den1=(z4temp(1:k*470)-mean(z1temp(1:k*470)))'*(z4temp(1:k*470)-mean(z1temp(1:k*470))); 
nrmse1(k)=sqrt(num1)/sqrt(den1);
end

for k=1:10
%Computer 4 with 2
num2=(z4temp(1:k*470)-z2temp(1:k*470))'*(z4temp(1:k*470)-z2temp(1:k*470)) ;
den2=(z4temp(1:k*470)-mean(z2temp(1:k*470)))'*(z4temp(1:k*470)-mean(z2temp(1:k*470))); 
nrmse2(k)=sqrt(num2)/sqrt(den2);
end

for k=1:10
%Computer 4 with 3
num3=(z4temp(1:k*470)-z3temp(1:k*470))'*(z4temp(1:k*470)-z3temp(1:k*470)); 
den3=(z4temp(1:k*470)-mean(z3temp(1:k*470)))'*(z4temp(1:k*470)-mean(z3temp(1:k*470))); 
nrmse3(k)=sqrt(num3)/sqrt(den3);
end

%==========================================================================
fprintf('%s\n','=====================================================');
fprintf('%s\n','-----------------------------------------------------');
fprintf('%s\n','NRMSE index result for different simulation intervals');
NRMSE=[NRMSE1; NRMSE2; NRMSE3; NRMSE4];
fprintf('%s\n',' Computer 1  Computer 2  Computer 3  Computer 4');
fprintf('%9.4f  %9.4f  %10.4f  %10.4f\n',NRMSE);
fprintf('%s\n','-----------------------------------------------------');
fprintf('%s\n','NRMSE index result assuming computer 4 as reference');
nrmse=[nrmse1; nrmse2; nrmse3];
fprintf('%s\n',' Computer 1  Computer 2  Computer 3 ');
fprintf('%9.4f  %9.4f  %10.4f\n',nrmse);
fprintf('%s\n','=====================================================');


%==========================================================================
%FIGURES
figure(1)
subplot(4,1,1),plot(t1,z1,'-','LineWidth',2,'Color',[1 0.6 0]);
xlabel('Time [s]','FontSize',16,'FontName','Times');
ylabel('Voltage [V]','FontSize',16,'FontName','Times');
set(gca,'fontsize',16,'FontName','Times') 
grid on;
box off;
subplot(4,1,2),plot(t2,z2,'-','LineWidth',2,'Color',[0 0.71 0.4]);
xlabel('Time [s]','FontSize',16,'FontName','Times');
ylabel('Voltage [V]','FontSize',16,'FontName','Times');
set(gca,'fontsize',16,'FontName','Times') 
grid on;
box off;
subplot(4,1,3),plot(t3,z3,'-','LineWidth',2,'Color',[0 0 1]);
xlabel('Time [s]','FontSize',16,'FontName','Times');
ylabel('Voltage [V]','FontSize',16,'FontName','Times');
set(gca,'fontsize',16,'FontName','Times') 
grid on;
box off;
subplot(4,1,4),plot(t4,z4,'-','LineWidth',2,'Color',[0 0 0]);
xlabel('Time [s]','FontSize',16,'FontName','Times');
ylabel('Voltage [V]','FontSize',16,'FontName','Times');
set(gca,'fontsize',16,'FontName','Times') 
grid on;
box off;
%--------------------------------------------------------------------------
figure(2)
plot(tC,zC,'-','LineWidth',2,'Color',[1 0 0]);
hold on
plot(t4,z4,'-','LineWidth',2,'Color',[0 0 0]);
xlabel('Time [s]','FontSize',20,'FontName','Times');
ylabel('Voltage [V]','FontSize',20,'FontName','Times');
set(gca,'fontsize',20,'FontName','Times')
grid on;
box off;
%--------------------------------------------------------------------------
figure(3)
plot(t1,z1,'LineWidth',5,'color',[1 0.6 0])
hold on
plot(t2,z2,'LineWidth',2,'color',[0 0.71 0.4    ])
hold on
plot(t3,z3,'LineWidth',2,'color',[0 0 1])
hold on
plot(t4,z4,'LineWidth',2,'color',[0 0 0])
xlabel('Time [s]','FontSize',20,'FontName','Times');
ylabel('Voltage [V]','FontSize',20,'FontName','Times');
set(gca,'fontsize',20,'FontName','Times') 
xlim([0.06 0.1])
grid on;
box off;
%--------------------------------------------------------------------------
figure(4)
plot(tC,zC,'-','LineWidth',2,'Color',[1 0 0]);
xlabel('Time [s]','FontSize',20,'FontName','Times');
ylabel('Voltage [V]','FontSize',20,'FontName','Times');
set(gca,'fontsize',20,'FontName','Times') 
grid on;
box off;