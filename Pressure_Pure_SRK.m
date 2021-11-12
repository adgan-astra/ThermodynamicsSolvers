% Calculating pressure for pure component using SRK
clear workspace;
clc;

R = 8.314;
T = 300 + 273.15;
vs = 5.226/10;
MW = 18;
vm = vs * MW / 1000;
Tc = 647.29;        % Kelvins
Pc = 22.09;         % MPa
Tr = T/Tc;  
w = 0.344;

a = 0.42747*R^2*Tc^2/(Pc*10^6);
b = 0.08664*R*Tc/(Pc*10^6);
alpha = (1+(0.48508+1.55171*w-0.15613*w^2)*(1-Tr^0.5))^2;

P = (R*T/(vm-b)*10^(-6)) - ((a*alpha/(vm*(vm+b)))*10^(-6));

fprintf (' P = %5.4f MPa \r\n', P);