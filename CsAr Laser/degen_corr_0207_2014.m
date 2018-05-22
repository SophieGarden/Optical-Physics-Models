%Modified by Wendy on Feb 7,2014. 
% Rate Equation Model to predict the dynamics of the CsAr Laser
% Predicts output with a guassian pump pulse
% Modified from v2 to be m^3/ns with a three-level model
% Paradigm shift -- just dump E_abs into system and see what pulse results
clear all;
clc;
clf;
%% Special Simulation Paramters
t_conv = 1e-10;  % factor to convert seconds into ns because minimum step size is stupidly fixed
magic_photons = 1; %13e-3; % this is added to the stimulated emission rate to keep the number of photons from going negative
T = 400;%K

%% General Physical constants
c = 2.99792458e8;  % m/s
h = 6.626068e-34;  % J-s
kb = 1.3806503e-23; % J/K
mp = 1.67262158e-27; % kg
q = 1.602e-19; % C

%% Cs constants
M_Cs = 133*mp;

g0 = 2;  % lower state property

gd1 = 2;
ld1 = 894.5929598610e-9; % wavelength in m % D1 line properties
td1 = 34.79190e-9; % lifetime

gd2 = 4;
ld2 = 852.3472758227e-9;  % wavelength in m % D2 line properties
td2 = 30.40577e-9; % lifetime
Ed2 = h*c/ld2; %energy of D2 photon
D_nu_D = sqrt(8*kb*T*log(2)/(M_Cs*c^2))*(c/ld2);
g_nu_0 = 2/D_nu_D*sqrt(log(2)/pi);
sigma_se_d2 = (1/td2)*ld2^2/(8*pi)*g_nu_0;
Is2 = Ed2/(sigma_se_d2*td2);

% t21 = 1e-10;  % lifetime of molecule

%% Empirical Cs data
a = dlmread('CsArKr.csv');  % absorption profile
wvl = linspace(830,850,100);
abs_spec = interp1(a(:,1),a(:,2),wvl);  % absorption coefficient k cm^-5
abs_spec = abs_spec*1e5; % absorption coefficient k m^-5 alpha = k[Cs][Ar]

% density vs temperature
P_Cs = (101325/760)*10^(2.881+4.165-3830/T);  % P_Cs is Cs pressure in Pa
N_Cs = P_Cs/kb*T;  % number density in m^-3


%% Experimental Parameters
% Physical Parameters
L  = 25*.0254; % .4;    % m   length of cavity
Lg = 4.125*.0254; % .2;   % m   length of gain medium
A  = (.003)^2; % area of pump
V  = L*A; % m^3 volume of laser cavity
Vg = Lg*A; % m^3 volume of gain medium intersecting with cavity mode

% Pump Parameters
tp = 1e-9; % pump FWHM
E_abs = 2; % pump energy %maybe mJ according to TG
sigma_tp = tp/sqrt(8*log(2));
pfctr = (E_abs/(Ed2*Vg))/(sigma_tp*sqrt(2*pi));  % MODIFY FOR PUMPING RATE 
% Precomputes exp prefactor so it doesn't have to be evaulated every time the function is called
expfctr = -(t_conv)^2/(2*sigma_tp^2);  % Precomputes exp argument so it " 
Rt = @(t) 1.5*pfctr*exp(expfctr*(t-50).^2);
pop_shift = E_abs/(Ed2*Vg);


% Cavity Parameters
% Q = 100;
% tp = Q*c/(2*pi*ld2);  % tau_p = Q/w0
R1 = .99;
R2 = .04;
tc = 2*L/(c*(1-R1*R2));
gamma_th = -1/(2*Lg)*log(R1*R2);
DN_th = gamma_th/sigma_se_d2;


%% Calculate Ouput Energy
degen = gd2/g0;
fc = 120;% gamma_0/gamma_th
N2_0 = (fc*DN_th+degen*N_Cs)/(1+degen);
R = 1.4*N2_0/td2;

phi = Is2*V/(c*Ed2);  % saturation of phi

% Factors are defined here to reduce multiplcations in function call
thefactor = sigma_se_d2*c/V;
thefactorV = thefactor*Vg;
A21 = 1/td2;
dN = @(t,N) ...
     t_conv*[+N(2)*A21+thefactor*N(3).*(N(2)-degen*N(1));... % -pfctr*exp(expfctr*t.^2)
      -N(2)*A21-thefactor*N(3).*(N(2)-degen*N(1));...
      1e6-N(3)/tc+thefactorV*(N(3)+magic_photons).*(N(2)-degen*N(1))];
% dN = @(t,N) ...
%      t_conv*[-R+N(:,2)*A21+thefactor*N(:,3).*(N(:,2)-degen*N(:,1));... % -pfctr*exp(expfctr*t.^2)
%       +R-N(:,2)*A21-thefactor*N(:,3).*(N(:,2)-degen*N(:,1));...
%       1-N(:,3)/tc+thefactorV*(N(:,3)+magic_photons).*(N(:,2)-degen*N(:,1))];

ode_opts = odeset('MaxStep',1e-11/t_conv);%,'Jacobian',JN); % 'NonNegative',[2]); % 'Vectorized','off'
% ode_opts = odeset('MaxStep',.01,'NonNegative',[3]); % 'Vectorized','off'
tstart = 0/t_conv;
tfin = 210e-9/t_conv;
[t,N] = ode23s(dN,[tstart tfin],[N_Cs-pop_shift pop_shift 1],ode_opts); % ode15s doesn't work

dN1 = -R+N(:,2)*A21+thefactor*N(:,3).*(N(:,2)-degen*N(:,1));
dN2 = R-N(:,2)*A21-thefactor*N(:,3).*(N(:,2)-degen*N(:,1));
dN3 = -N(:,3)/tc+thefactorV*(N(:,3)+magic_photons).*(N(:,2)-degen*N(:,1));

%% Plot Results
Inversion = N(:,2)-degen*N(:,1);
Output_Power = N(:,3)*Ed2/(tc);  % units of photons/ns

% Plot Inversion
figure;
semilogy(t,(Inversion),':')
hold on;
semilogy(t,abs(N(:,2)),'r')
semilogy(t,abs(N(:,1)),'g')
line([t(1) t(end)],[DN_th DN_th],'LineStyle','--');
line([t(1) t(end)],[N_Cs N_Cs],'LineStyle','--');
%hold off;

% Plot Output Power
figure;
semilogy(t,abs(N(:,3)))
% ylim([0 1.1*max(N(:,2))]);
title('N_p');

E_in = Ed2*Vg*trapz(t_conv*t,Rt(t));
E_out = trapz(t_conv*t,Output_Power);
a = thefactor*N(:,3).*(N(:,2)-degen*N(:,1));
photon_lifetime = 1e9./(1/tc-thefactorV*(N(:,2)-degen*N(:,1)));
Photons_Out = trapz(t_conv*t,a);
fprintf('Input Energy E = %2.2e J\n',E_in)
fprintf('Output Energy E = %2.2e J\n',E_out)
fprintf('Net Photons created = %2.2e\n',Photons_Out)
title('Energy E ');
%% CW Output Power
P_out =  A*Is2*(1-R2)*(Lg*1.5*gamma_th*fc+log(R1*R2))./((1-sqrt(R1*R2)).*(1+sqrt(R2/(R1))));
fprintf('Rigrod Power P = %2.2f\n',P_out)
fprintf('Simulated Power P = %2.2f\n',Output_Power(end))
title('CW Output Power');
% % % %% Extra calculations 
% % % % Theoretical SS Output Power
% % % DN_0 = 4/3*extra*N_Cs;
% % % Pout = A*Is2/2*(DN_0/DN_th-1);
% % % Photons_SS = Pout*tc/Ed2;
% % % % rabi frequency
% % % % tc - mean time between collisions
