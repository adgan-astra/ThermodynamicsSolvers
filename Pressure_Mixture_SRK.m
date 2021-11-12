
% Calculating pressure for mixtures using SRK
clear workspace;
clc;

y1 = 0.5;
y2 = 0.5;
MW1 = 44.01;
MW2 = 58.12;
MW = y1*MW1+y2*MW2;

R = 8.314;
T = 377.5;
vs = 14.68*10^(-3)/1.0;      % m3/kg
vm = vs * MW / 1000;         % m3/gmol
kij = 0.5; 

w1 = 0.225;
w2 = 0.193;

Tc1 = 304.19;             % Kelvins
Pc1 = 7.382;         % MPa
Tr1 = T/Tc1; 

Tc2 = 425.16;             % Kelvins
Pc2 = 3.796;         % MPa
Tr2 = T/Tc2; 

a1 = 0.42747*R^2*Tc1^2/(Pc1*10^6);
b1 = 0.08664*R*Tc1/(Pc1*10^6);
alpha1 = (1+(0.48508+1.55171*w1-0.15613*w1^2)*(1-Tr1^0.5))^2;
aam1 = a1*alpha1;

a2 = 0.42747*R^2*Tc2^2/(Pc2*10^6);
b2 = 0.08664*R*Tc2/(Pc2*10^6);
alpha2 = (1+(0.48508+1.55171*w2-0.15613*w2^2)*(1-Tr2^0.5))^2;
aam2 = a2*alpha2;

aam = y1*y1*aam1 + y2*y2*aam2 + 2*y1*y2*(1-kij)*sqrt(aam1*aam2);
bm = y1*b1+y2*b2;

P = (R*T/(vm-bm)*10^(-6)) - ((aam/(vm*(vm+bm)))*10^(-6));
z = (P*10^6)*vm/(R*T);             % Compressibility Factor

fprintf (' \n Pressure and Compressibility factor from RK EOS: \r');
fprintf (' P = %8.5f MPa \r', P);
fprintf (' z = %8.5f  \r\n', z);
