 %Rayleigh scattering
pia = 3.14;
d = 100e-9; %diameter of quantum dot
lamd = 532e-9; %pump wavelength
n1 = 1.4973; %solution Tolluene index
n2 = 2.59748; %quantum dot index
n = n2/n1;%relative index
m = 5e-6; % Kg
V = 1e-6;%m^3
Na = 6.02214129e23;%Avogadro's number)
Am = 112.4+78.96; %Cd + Se
sgmS = 2*pia^5/3*d^6/lamd^4*((n^2-1)/(n^2+2))^2;%scattering cross-section

N = m/Am*Na/V;%concentration

L0 = 1/(sgmS*N); % in SI base unit meter
L0inMicron = L0*1e6 % in unit Micron meter