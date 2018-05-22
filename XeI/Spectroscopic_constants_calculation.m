% spectraconstants

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

%-----------------------------

makeparams();
load('SimParams');
T=300; % kelvin 
%----------------For XeI----------------------
M_Xe = 131.293; %u
M_I = 126.90447;% u
M_red = 1/(1/M_Xe+1/M_I);%
%-------------------------------------
rmin  =  2.7;% %Angstrom
rstep =  5e-4*1;
rmax  =  4.7;%
R = (rmin:rstep:rmax)';%prime denotes the transpose of a function %R is a single column vector

%---U_B-----
U_B = @(R) 73160.89 + b_B*exp(-beta_B*R)-1.16141e5*R.^(-1)-3.9e4*R.^(-3)-5.3e5*R.^(-4)-6.3e5*R.^(-6);
numbound = 60;

[wvfns_u, energies_u] = Bound_Wvfns_Mat_WendyV5(R,U_B(R),M_red,pc,numbound); %wavfns' columns corresponding eigenvectors.

[Te,Imin]=min(U_B(R));
De=U_B(1e8)-Te;
re=R(Imin);
Te_sensitivity=max(U_B(re-rstep)-Te,U_B(re-rstep)-Te);
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[200,0.5],...
               'StartPoint',[105 0.2]);
ft = fittype('we*(v+0.5)-wexe*(v+0.5)^2','independent','v','coefficients',{'we','wexe'},'options',fo);

v=(0:1:numbound-1)';
G=energies_u-Te;
[curve2,gof2] = fit(v,G,ft)
%%
figure;
plot(curve2,v, G);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---U_X-----
U_L_Morse = @(R) De_l*((exp(-beta_l*(R-Re_l))-1).^2-1);
U_R_Morse = @(R) De_r*((exp(-beta_r*(R-Re_r))-1).^2-1); 

U_X = @(R) U_L_Morse(R)./(1+exp(a_off_l*(R-Rs_off_l)))...
    +U_R_Morse(R)./((1+exp(-a_on_l*(R-Rs_on_l))).*(1+exp(a_off_r*(R-Rs_off_r))))...
-(1.73e6*R.^(-6)+1.21e7*R.^(-8))./(1+exp(-a_on_r*(R-Rs_on_r))); 

[Te_x,Imin_x]=min(U_X(R));
De_x=0-Te_x

re_x=R(Imin_x)
DeX_sensitivity=max(abs(U_X(re_x-rstep)-Te_x),abs(U_X(re_x-rstep)-Te_x))
UXRe=U_X(re)
UXRe=(U_X(re+rstep)-U_X(re-rstep))/(2*rstep)
