% Calculating pressure for pure component using RK
clear workspace;
clc;

R = 8.314;
T = 250 + 273.15;
vs = 24.06/10;  % m3/kg
MW = 18;
vm = vs * MW / 1000;  % m3/gmol
Tc = 647.29;        % Kelvins
Pc = 22.09;         % MPa

a = 0.42748*R^2*Tc^2.5/(Pc*10^6);
b = 0.08664*R*Tc/(Pc*10^6);

P = (R*T/(vm-b)*10^(-6)) - ((a/(T^0.5*vm*(vm+b)))*10^(-6));

fprintf (' P = %5.4f MPa \r\n', P);