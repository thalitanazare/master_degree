% x=linspace(0,10,1000);
% y=2*sin(2*pi*60*x);

load Jerk_1cond.txt 
load Jerk_2cond.txt 
load Jerk_3cond.txt 
load Jerk_4cond.txt 
load Jerk_circ.txt %voltage data collected in the real circuit.

load Jerk4.txt
load Jerk17.txt


z1 = Jerk_1cond(:,2);
z2 = Jerk_2cond(:,2); 
z3 = Jerk_3cond(:,2);
z4 = Jerk_4cond(:,2);
zC = Jerk_circ(:,1);

zLT4 = Jerk4(:,1);
zLT17 = Jerk17(:,1);

%Convertendo para uma escala de 0 a 255
y1=round(255*(z1-min(z1))./(max(z1)-min(z1)));
y2=round(255*(z2-min(z2))./(max(z2)-min(z2)));
y3=round(255*(z3-min(z3))./(max(z3)-min(z3)));
y4=round(255*(z4-min(z4))./(max(z4)-min(z4)));
yC=round(255*(zC-min(zC))./(max(zC)-min(zC)));

yLT4=round(255*(zLT4-min(zLT4))./(max(zLT4)-min(zLT4)));
yLT17=round(255*(zLT17-min(zLT17))./(max(zLT17)-min(zLT17)));

E1=myentropy(y1(:))
E2=myentropy(y2(:))
E3=myentropy(y3(:))
E4=myentropy(y4(:))
EC=myentropy(yC(:))

ELT4=myentropy(yLT4(:))
ELT17=myentropy(yLT17(:))