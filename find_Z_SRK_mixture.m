% SRK EOS - Mixture
% Newton-Raphson method to solve for Z using SRK
clear workspace; clc;

y1 = 0.5;
y2 = 0.5;
MW1 = 30.07;
MW2 = 42.08;
MW = y1*MW1+y2*MW2;

R = 8.314;
T = 277.59;   % Kelvins
P = 1;        % MPa
kij = 0.05; 

w1 = 0.22394;
w2 = 0.15238;

Tc1 = 305.4;                   % Kelvins
Pc1 = 4883865*10^(-6);         % MPa
Tr1 = T/Tc1; 

Tc2 = 365.0;                   % Kelvins
Pc2 = 4620420*10^(-6);         % MPa
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

Am = aam*P*(10^6)/(R^2*T^2);
Bm = bm*P*(10^6)/(R*T);

% Coefficients
a = 1.0;                   % for Z^3
b = -1.0;                  % for Z^2 
c = Am-Bm-(Bm*Bm);             % for Z
d = -Am*Bm;%                 % constant term

TOL = 1;   % Initial large tolerance for while loop
Zold = 1;  % Initial Guess
i = 0;     % Counter

fprintf('\r\n');
% Print formatting
fprintf('%5s %10s %10s %14s %8s\r\n', 'i', 'Z2', 'F(Z2)', 'Fprime(Z2)','Error');
fprintf('%5d  %9.2f   %10.4f  %10.4f \t\t-- \r\n',i, Zold,fun(Zold,a,b,c,d),fprime(Zold,a,b,c));

while TOL > 1e-6  % run until tolerance is less than 1e-6
    Znew = Zold - fun(Zold,a,b,c,d)/fprime(Zold,a,b,c);
    TOL  = abs(Znew-Zold);
    error = abs(Znew-Zold);
    i = i + 1;
    fprintf('%5d  %9.6f  %11.6f  %10.4f %10.5f \r\n',i, Znew,fun(Znew,a,b,c,d),fprime(Znew,a,b,c), error);
    Zold = Znew;
end

% Molar volume
vmolar = Znew*R*T/(P*10^6);
v = vmolar*1000/MW;
fprintf('  molar v = %8.5f m3/gmol\r\n',vmolar);
fprintf('  specific v = %8.5f m3/kg \r\n',v);

% Polynomial form in Z
function [f_t] = fun(Z,a,b,c,d)
    f_t = a*Z^3 + b*Z^2 + c*Z + d;
end

% Z prime
function [fp] = fprime(Z,a,b,c)
    fp = 3*a*Z^2 + 2*b*Z + c;
end
