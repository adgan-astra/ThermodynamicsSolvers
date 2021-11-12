% Calculating pressure for mixtures using RK
% The following is the solution to problem 1.2 and 1.3 
% in the Phi_Phi_Problem notes

clear workspace;
clc;

y1 = 0.65;
y2 = 0.35;
MW1 = 30.07;
MW2 = 42.08;
MW = y1*MW1+y2*MW2;            % Average Molecular Weight

R = 8.314;                     % m3*Pa/gmol*K
T = 280.0;
vs = 27.5*10^(-3)/0.5;         % m3/kg
vm = vs * MW / 1000;           % m3/gmol

Tc1 = 305.4;                   % Kelvins
Pc1 = 4.883865;                % MPa

Tc2 = 365.0;                   % Kelvins
Pc2 = 4.620420;                % MPa

a1 = 0.42747*R^2*Tc1^2.5/(Pc1*10^6);
b1 = 0.08664*R*Tc1/(Pc1*10^6);

a2 = 0.42747*R^2*Tc2^2.5/(Pc2*10^6);
b2 = 0.08664*R*Tc2/(Pc2*10^6);

am = y1*y1*a1 + y2*y2*a2 + 2*y1*y2*sqrt(a1*a2);
bm = y1*b1+y2*b2;
P = (R*T/(vm-bm)*10^(-6)) - ((am/(T^0.5*vm*(vm+bm)))*10^(-6)); % MPa
z = (P*10^6)*vm/(R*T);         % Compressibility Factor

fprintf (' \n Pressure and Compressibility factor from RK EOS: \r');
fprintf (' P = %8.5f MPa \r', P);
fprintf (' z = %8.5f  \r\n', z);

% Now find fugacity coefficient
B1 = b1*(P*10^6)/(R*T);
B2 = b2*(P*10^6)/(R*T);
Bm = bm*(P*10^6)/(R*T);

A1 = a1*(P*10^6)/(R^2*T^2.5);
A2 = a2*(P*10^6)/(R^2*T^2.5);
Am = am*(P*10^6)/(R^2*T^2.5);

LN_Phi1v = B1/Bm*(z-1)-log(z-Bm)+ (Am/Bm*(B1/Bm-(2*sqrt(A1/Am)))*log(1+Bm/z));
Phi1v = exp(LN_Phi1v);

LN_Phi2v = B2/Bm*(z-1)-log(z-Bm)+ (Am/Bm*(B2/Bm-(2*sqrt(A2/Am)))*log(1+Bm/z));
Phi2v = exp(LN_Phi2v);

fprintf (' Fugacity coefficient of component 1: \r');
fprintf (' ln Phi_1_v = %7.5f  \r', LN_Phi1v );
fprintf (' Phi_1_v = %7.5f  \r\n', Phi1v );

fprintf (' Fugacity coefficient of component 2: \r');
fprintf (' ln Phi_2_v = %7.5f  \r', LN_Phi2v );
fprintf (' Phi_2_v = %7.5f  \r\n', Phi2v );

