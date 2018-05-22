clear;
clc

m=200.59*1.66053892e-27;%kg 
m1=m*2;%Hg dimer
m2=m*3;% Hg trimer
c=3e8;% [m/s]
lamda1 = 335e-9; % [m]
lamda2 = 485e-9; % [m]
v01=c/lamda1; %[Hz]
v02=c/lamda2; %[Hz]
k = 1.3806488e-23; %Blotzamann constant in unit of  m2 kg s-2 K-1
T=500;%[K]
kT = k*T; %[m2 kg s-2]
deltaVd1 = 2*v01/c*sqrt(2*log(2)*kT/m1);%[Hz]
deltaLamda1 = c/v01^2*deltaVd1; %[meter]
deltaVd2 = 2*v02/c*sqrt(2*log(2)*kT/m2);%[Hz]
deltaLamda2 = c/v02^2*deltaVd2; %[meter]
A335 = 1e6;%[sec-1]

% guess the population inversion
deltaN = 1e17*1e-6; %[m^-3]

V = 3.14*2.54^2/4*4*1e-6;
n = 1; %refractive index
gv01 = 2/deltaVd1*sqrt(log(2)/3.14);
sigma1 = A335*lamda1^2/(8*3.14*n^2)*gv01; %cross section
r1 = sigma1*deltaN*V;%gain
gv02 = 2/deltaVd2*sqrt(log(2)/3.14);
sigma2 = A335*lamda1^2/(8*3.14*n^2)*gv02;
r2 = sigma2*deltaN*V;
