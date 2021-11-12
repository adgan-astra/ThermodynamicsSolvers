% SRK EOS - Pure Component
% Newton-Raphson method to solve for Z using SRK
clear workspace;

TOL = 1;   % Initial large tolerance for while loop
Zold = 1;  % Initial Guess
i = 0;     % Counter

% Input
T_C = 350;           % Temperature in Celcius
T = T_C + 273.15;    % Kelvins
P = 3;               % MPa
w = 0.344;           % Accentric factor

R = 8.314;
Tc = 647.29;  % Kelvins
Pc = 22.09;   % MPa
Tr = T/Tc;   

a = 0.42747*R^2*Tc^2/(Pc*10^6);
b = 0.08664*R*Tc/(Pc*10^6);
alpha = (1+(0.48508+1.55171*w-0.15613*w^2)*(1-Tr^0.5))^2;

A = a*alpha*P*(10^6)/(R^2*T^2);
B = b*P*(10^6)/(R*T);

% Coefficients
a = 1.0;                   % for Z^3
b = -1.0;                  % for Z^2 
c = A-B-(B*B);             % for Z
d = -A*B;%                 % constant term

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
MW = 18;
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