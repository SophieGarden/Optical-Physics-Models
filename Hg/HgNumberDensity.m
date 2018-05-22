% Hg number density-temperature 
%Mercury vapor pressure
%
clear all;
clc;

NA = 6.02214129e23;% mol?1, Avogadro constant 
R = 8.314472; % J/(mol K)
%m = 200.59; %g/mol

a1 = -4.57618368;%0.0472
a2 = -1.40726277; %0.8448
a3 = 2.36263541; %0.8204
a4 = -31.0889985; %1.3439
a5 = 58.0183959; %2.4999
a6 = -27.6304546; %1.1798
Tc = 1764; %K
Pc = 167; %Mpa

i = 1;
for n = [0.1,0.4,0.7,1,1.5,2,3,4,6,8,10,20,30,40,50,60,70,80,90,100,200]*1e17*1e6      %m^-3

rho = n/NA;
syms T
P = rho*R*T*1e-6;%Mpa
tau = 1-T/Tc;

solT = vpasolve(log(P/Pc) == (Tc/T)*(a1*tau+a2*tau^1.89+a3*tau^2+a4*tau^8+a5*tau^8.5+a6*tau^9),T,[-Inf Inf]);
celsiusT(1,i) = double(solT) - 273.15;
celsiusT(2,i) = n/1e6;
nn(i) = n;
i = i+1;
end
