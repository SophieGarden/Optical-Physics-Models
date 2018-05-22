clear;
makeparams();
load('SimParams');

% % % % % July 12th DATA. wendy version 


%===Upper B state===
b_B = 1.0183e7;
beta_B = 2.160;

%===Ground X state===
% Left morse,inner wall:
De_l=242;
beta_l=0.850;
Re_l=4.504;

% Switch functions, joint
a_off_l=8.1;
Rs_off_l=3.562;
a_on_l=3.65;
Rs_on_l=3.56;

%Right morse, well
De_r=172;
beta_r=0.9;
Re_r=4.02;

%Switch functions, Right slopes
a_off_r=2.08;
Rs_off_r=5.34;
a_on_r=1.46;
Rs_on_r=6.4;


%---------------------
T=300; % kelvin 
rmin  =  2.75;%-0.5; %Angstrom
rstep =  5e-4*1;
rmax  =  10;%+0.5;

%----------------------------------------------
R = (rmin:rstep:rmax)';%prime denotes the transpose of a function %R is a single column vector
R_far = (rmin:rstep:rmax+300-200-50)';% R_far is for calculation of the value of the wavefunction at the very far end
% Iterate through angular moentum quantum numbers


% Ground State Energies
% Energy stepsize
E_step = 4;  % cm^-1

% Upper State Energies
numbound = 60;
%----------------For XeI----------------------
M_Xe = 131.293; %u
M_I = 126.90447;% u

%----------------For KrF----------------------
% M_Xe = 83.798; %u
% M_I = 18.998403;% u

%----------------M_red----------------------
M_red = 1/(1/M_Xe+1/M_I);

U_B = @(R) 73160.89 + b_B*exp(-beta_B*R)-1.16141e5*R.^(-1)-3.9e4*R.^(-3)-5.3e5*R.^(-4)-6.3e5*R.^(-6);
%---U_X-----
U_L_Morse = @(R) De_l*((exp(-beta_l*(R-Re_l))-1).^2-1);
U_R_Morse = @(R) De_r*((exp(-beta_r*(R-Re_r))-1).^2-1); 

U_X = @(R) U_L_Morse(R)./(1+exp(a_off_l*(R-Rs_off_l)))...
    +U_R_Morse(R)./((1+exp(-a_on_l*(R-Rs_on_l))).*(1+exp(a_off_r*(R-Rs_off_r))))...
-(1.73e6*R.^(-6)+1.21e7*R.^(-8))./(1+exp(-a_on_r*(R-Rs_on_r))); 

potJ = @(R) (R*1e-10).^(-2)*(pc.hbar)^2/(2*M_red*pc.mp*pc.icm_2_J);%1e-10 is for angstrom unit

%rotational quantum number
J_min  = 0;
J_step = 25;%5;
J_max  =  325;
j=1;
for J = J_min:J_step:J_max  % iterate through rotational quantum numbers
    % add in effects of rotation to potential
    potpts_u(:,j) = U_B(R) + J*(J+1)*potJ(R);
    potpts_g(:,j) = U_X(R) + J*(J+1)*potJ(R);
    potpts_g_far(:,j) = U_X(R_far) + J*(J+1)*potJ(R_far);
j=j+1;
end
csvwrite('potpts_u.dat',horzcat(R,potpts_u));
csvwrite('potpts_g.dat',horzcat(R,potpts_g));

%%
figure;
plot(R, potpts_u);hold on;

xlabel('Interatomic Distance (Angstrom)');
ylabel('Energy (cm^{-1})');
% legend('');
title('B State Potential Curve of XeI');
xlim([2.7 4.7]);
%%
figure;
plot(R_far, potpts_g_far);hold on;

xlabel('Interatomic Distance (Angstrom)');
ylabel('Energy (cm^{-1})');
% legend('');
title('X State Potential Curve of XeI');
xlim([2.7 100]);