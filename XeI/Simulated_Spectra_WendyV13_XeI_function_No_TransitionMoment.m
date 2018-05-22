% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Wendy Chen LOPE
% % V 13 05/28/2017 Add all the 16 parameters to the function
% %  V12 07/22/2016 Fine tune U_X and recalculate Transiiton Moment
% %   V11 07/20/2016, Add a condition for  setting the upper limit of J,ie
% % J_max
% % V10 07/15/2016, Correct the codes for setting the lower energy limit
% % for free states U_X. The peaks match much better. The peak positions at
% % larger V' is determined mostly by U_X, instead of U_B.
% % V9 4/14/2016, for a certain set of J_min:J_step:J_max, run Matrix Numerov method first, record the values E_u_MN,
% %  Run Numervo-Cooley method, record E_u, then E_u could be the initial
% %  values of energies. Not sure how much difference the two methods can
% % have though. Numerove-Cooley methods generates a more smoother
% % wavefunction.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function FB_Absportion(T)
% T - temperature in Kelvin
% This function returns an simulated spectra

function [EEE,SSS] = Simulated_Spectra_WendyV13_XeI_function_No_TransitionMoment(b_B,beta_B,De_l,beta_l,Re_l,De_r,beta_r,Re_r,a_off_l,Rs_off_l,a_on_l,Rs_on_l,a_off_r,Rs_off_r,a_on_r,Rs_on_r)

%% Parameters for Simulation
% Load in the parameters of the elements

makeparams();
load('SimParams');
T=300; % kelvin 
% Range of distance between molecules and number of points
%----------------For XeI----------------------
rmin  =  2.7;%-0.5; %Angstrom
rstep =  5e-3*1;
rmax  =  4.7;%+0.5;

%----------------------------------------------
R = (rmin:rstep:rmax)';%prime denotes the transpose of a function %R is a single column vector
R_far = (rmin:rstep:rmax+300-200-50)';% R_far is for calculation of the value of the wavefunction at the very far end
% Iterate through angular moentum quantum numbers
%rotational quantum number
J_min  = 0;
J_step = 5;%5;
J_max  =  320*1;

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
M_red = 1/(1/M_Xe+1/M_I);%

% dd = pc.h./sqrt((8*pc.k*T)*2*pc.mp*M_red)*1e10; %max grid spacing 0.0783 A

%% Load in Potential Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % 16 parameters in total, [2 for U_B: b_B, beta_B], [14 for U_X:
% De_l,beta_l,Re_l,De_r,beta_r,Re_r,a_off_l,Rs_off_l,a_on_l,Rs_on_l,a_off_r,Rs_off_r,a_on_r,Rs_on_r]
% % b_B = 1.0169e7;
% % beta_B=2.16;
% % De_l=240;
% % beta_l=0.852;
% % Re_l=4.50;
% % De_r=182;
% % beta_r=0.870;
% % Re_r=3.966;
% % a_off_l=8;
% % Rs_off_l=3.57;
% % a_on_l=3.8;
% % Rs_on_l=3.6;
% % a_off_r=2;
% % Rs_off_r=5.3;
% % a_on_r=1.5;
% % Rs_on_r=6.5;
%---U_B-----
U_B = @(R) 73160.89 + b_B*exp(-beta_B*R)-1.16141e5*R.^(-1)-3.9e4*R.^(-3)-5.3e5*R.^(-4)-6.3e5*R.^(-6);
%---U_X-----
U_L_Morse = @(R) De_l*((exp(-beta_l*(R-Re_l))-1).^2-1);
U_R_Morse = @(R) De_r*((exp(-beta_r*(R-Re_r))-1).^2-1); 

U_X = @(R) U_L_Morse(R)./(1+exp(a_off_l*(R-Rs_off_l)))...
    +U_R_Morse(R)./((1+exp(-a_on_l*(R-Rs_on_l))).*(1+exp(a_off_r*(R-Rs_off_r))))...
-(1.73e6*R.^(-6)+1.21e7*R.^(-8))./(1+exp(-a_on_r*(R-Rs_on_r))); 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % % % % % XeI FANG SHEN's Thesis pg 76
% % % % % U_B = @(R) 73160.89 + 1.0169e7*exp(-2.16*R)-1.16141e5*R.^(-1)-3.9e4*R.^(-3)-5.3e5*R.^(-4)-6.3e5*R.^(-6);%Shen's
% % % % % 
% % % % % %%%%%%%%%  adjusted   %%%%%%%%%
% % % % % U_L_Morse = @(R) 240*((exp(-0.852*(R-4.500))-1).^2-1);%%Shen's
% % % % % U_R_Morse = @(R) 182*((exp(-0.870*(R-3.966))-1).^2-1);  
% % % % % U_X = @(R) U_L_Morse(R)./(1+exp(8*(R-3.57)))...
% % % % %     +U_R_Morse(R)./((1+exp(-3.8*(R-3.6))).*(1+exp(2*(R-5.3))))...
% % % % % -(1.73e6*R.^(-6)+1.21e7*R.^(-8))./(1+exp(-1.5*(R-6.5))); %Shen's

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XeI potential curve_RB Jones potential curve
% U_B =@(R) 73130 + 8.4346e6*exp(-2.145741*R)-1.16141e5*R.^(-1)-5.3e5*R.^(-4);
% U_X =@(R) 3.5e7*R.^(-9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
potJ = @(R) (R*1e-10).^(-2)*(pc.hbar)^2/(2*M_red*pc.mp*pc.icm_2_J);%1e-10 is for angstrom unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transition moment Fang Shen's Thesis pg 81
% mur = @(R) 1.25^2./((1.25^2+(R-3.33).^2).*(1+exp(20*(R-4.0))))+ 0.82^2./((0.82^2+(R-3.40).^2).*(1+exp(-20*(R-4.0))));
% murpts = mur(R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Adjusted Transition moment for XeI
% k = 0.82*1.5;
% m = 0.4;
% s = 0.5;
% ss=0.2;
% kk=20;% 20
% mur = @(R) (1.25*s)^2./((1.25*s^2+(R-3.33-m).^2).*(1+exp(kk*(R-4.0-ss))))+ (k)^2./((k^2+(R-3.40).^2).*(1+exp(-kk*(R-4.0-ss))));
% murpts = mur(R);
% murpts = murpts ./max(murpts);
%%%%%%%%%%%%%%%%%%%%%
%directly calculated from equation MuV=(Hexp./Hsim).^(.5);
% mur = @(R) .4545^2./((.4545^2+(R-3.703).^2).*(1+exp(20*(R-4.00))))+.4187^2./((.4187^2+(R-3.771).^2).*(1+exp(-20*(R-4.0))));

% 
%        a1 =      0.4478 ;% (0.4209, 0.4747)
%        b1 =       3.701 ;% (3.69, 3.712)
%        a2 =      0.3711 ;% (0.324, 0.4182)
%        b2 =        3.81  ;%(3.764, 3.856)
% mur = @(R) a1^2./((a1^2+(R-b1).^2).*(1+exp(20*(R-4.0))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-20*(R-4.0))));
%%%%%
%        a1 =      0.4666 ;% (0.4488, 0.4844)
%        b1 =       3.707 ;%  (3.699, 3.715)
%        c =       4.198 ;% (4.132, 4.263)
%        a2 =      0.1929 ;% (0.1299, 0.2558)
%        b2 =       4.022 ;% (3.944, 4.1)
% mur = @(R) a1^2./((a1^2+(R-b1).^2).*(1+exp(20*(R-c))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-20*(R-c))));
% murpts = mur(R);
%%%%%%%%%%%%%%%%%%%%%%%
% % % % %% calculated from equation MuV=(Hexp./Hsim).^(.5);    
% % % % %        a1 =      0.4679 ;% (0.4493, 0.4866)
% % % % %        b1 =       3.707 ;% (3.699, 3.715)
% % % % %        d =       24.99  ;%(0.2788, 49.69)
% % % % %        c =       4.176  ;%(4.085, 4.268)
% % % % %        a2 =       0.208 ;% (0.1288, 0.2873)
% % % % %        b2 =       4.003 ;% (3.903, 4.102)
% % % % %  %  adjusted
% % % % %        a1 =      0.4679 +0.05;% (0.4493, 0.4866)
% % % % % b1 =       3.707 +0.01;% (3.699, 3.715)
% % % % % d =       24.99  ;%(0.2788, 49.69)
% % % % % c =       4.176  ;%(4.085, 4.268)
% % % % % a2 =       0.208+0.05 ;% (0.1288, 0.2873)
% % % % % b2 =       4.003 +0.01;% (3.903, 4.102)
% % % %        a1 =      0.5737 ;%  (0.5513, 0.5961)
% % % %        b1 =       3.704  ;% (3.692, 3.715)
% % % %        d =       41.57 ;%  (-36.61, 119.7)
% % % %        c =       4.282 ;%  (4.083, 4.482)
% % % %        a2 =      0.2771;%   (-0.3736, 0.9278)
% % % %        b2 =       3.948  ;% (2.939, 4.958)
% % % %  mur = @(R)  a1^2./((a1^2+(R-b1).^2).*(1+exp(d*(R-c))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-d*(R-c))));
% % % % murpts = mur(R);
%%%%%%%%%%%%%%%%%%%%%
% %  constant mur
% murpts = 1;
 murpts = ones(length(R),1);
%%%%%%%%%%%%%%%%%%%%%
%% Set up the x (wavelength) axis
nnnn=(J_max-J_min)/J_step+1;
potpts_u = zeros(length(R),nnnn);
potpts_g = zeros(length(R),nnnn);
potpts_g_far = zeros(length(R_far),nnnn);
E_min = zeros(1,nnnn);
E_max = zeros(1,nnnn);
E_u = zeros(numbound,nnnn);
W_u = zeros(length(R), numbound,nnnn);
E_u_MN = zeros(numbound,nnnn);
W_u_MN = zeros(length(R), numbound,nnnn);
DeB = zeros(1,nnnn);
% DeBR =zeros(1,nnnn);
TeB1 = zeros(1,nnnn);
% TeB2 = zeros(1,nnnn);
% TeB1R = zeros(1,nnnn);
% TeB2R =zeros(1,nnnn);
Rwell=zeros(1,nnnn);
locs2=zeros(1,nnnn);

DIF=zeros(length(R_far)-1,nnnn);

EE0 = zeros(numbound*462,nnnn);
SS0 = zeros(numbound*462,nnnn);
%% Simulate the spectrum
j=1;
for J = J_min:J_step:J_max  % iterate through rotational quantum numbers
    % add in effects of rotation to potential
    potpts_u(:,j) = U_B(R) + J*(J+1)*potJ(R);
    potpts_g(:,j) = U_X(R) + J*(J+1)*potJ(R);
    potpts_g_far(:,j) = U_X(R_far) + J*(J+1)*potJ(R_far);
    %+++++++++++++Numerov-cooley+++++++++++++++++++++++++++++++++++++++++++++++++++
    %%%%%---------------------Numerov-cooley
    %%%%%%%----------version included dependence of J for initial energy
    method =2;
    if method ==1;% indicator for Numerov-cooley method
        load('E_u_MN.mat');
        %load('E_u.mat'); %For second time run N-C method
        energies0 = E_u(:,j);% !J_min, J_max,J_step has to be the same with the one generates E_u
        %                   load('E_u_MN.mat');
        %                energies0 = E_u_MN(:,j);% !J_min, J_max,J_step has to be the same with the one generates E_u
        
        energies0 = sort(energies0); % energies0 in ascending order, lower to higher energies
        [wvfns_u, energies_u] = Bound_Wvfns_Wendy_NCV9(R,potpts_u(:,j),M_red,pc,numbound,energies0); %wavfns whose columns are the corresponding eigenvectors.
        E_u(:,j) = energies_u; %record values, for future assignment of initiate values of energies0
        W_u(:,:,j) = wvfns_u;
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %%%%%%---------- Matrix numerov------------------------------------------------
    elseif method == 2; % indicator for Matrix numerov method
        [wvfns_u, energies_u] = Bound_Wvfns_Mat_WendyV5(R,potpts_u(:,j),M_red,pc,numbound); %wavfns' columns corresponding eigenvectors.
        E_u_MN(:,j) = energies_u;
        W_u_MN(:,:,j) = wvfns_u;
    end
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % set the lower energy limit for the free ground state
    % thermally free energy is ~ 300 cm^-1
    %     if min(potpts_g(:,j))<0
    %         E_min(j) = 1+3; % E_min = 0 wont work in the code
    %     else
    %         E_min(j) = min(potpts_g(:,j))+1+3; % +1 is for fun
    %     end
    
    % % % free collision paris should have energy larger than De'' ( Dissociation
    % % % energy of lower state)
    %     [TeB1(j),TeB1R(j)]= min(potpts_g(:,j)); % miminum of the interested range
    %     [DeB(j),DeBR(j)]= max(potpts_g_far(R>4.0,j)); % could be De of B state, or just local maximum if R>4.05 is monotonic decaying
    %     [TeB2(j),TeB2R(j)]= min(potpts_g(1:DeBR(j),j));% local minimum
    %     if  TeB2(j) < DeB(j) && DeB(j)>0
    %         E_min(j) = DeB(j)+1+3; % E_min = 0 wont work in the code
    %     elseif TeB2(j) < DeB(j) && DeB(j)<=0
    %         E_min(j) = 0 +1+3;
    %     elseif TeB2(j) >= DeB(j)
    %         E_min(j) = TeB1(j)+1+3; % +1 is for fun
    %     end
    
    % % % free collision paris should have energy larger than De'' ( Dissociation  energy of lower state)
    [TeB1(j),~]= min(potpts_g(:,j)); % miminum of the interested(calculated) range
    DIF(:,j) = diff(potpts_g_far(:,j));
    R5=find(R_far==5);%according to results seen from potential_curve.m, if R>5, the U_X is monotonic decay
    [~,locs2(j)] = min(abs(DIF(1:R5,j))); %locs2 is the indices of R of well
    % Position1(i)=locs2;
    Rwell(j)=R_far(locs2(j));
    DeB(j)= max(potpts_g_far(locs2(j):end,j));%barrier height
    %Emin=0 won't work for this code
    if Rwell(j)<5 && DeB(j)>0
        E_min(j) = DeB(j)+0;%%%increase barrier height %35 seems good. % 30,40 works;10, 20,50,60,80,100 doesn't work
    elseif Rwell(j)<5 && DeB(j)<=0
        E_min(j) =0;
    elseif Rwell(j)>=5
        E_min(j) =TeB1(j);
    else
        msg = 'Error occurred.';
        error(msg)
    end
    if E_min(j)<10
        E_min(j)=10;
    end
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % Free state wavefunctions calculation
    f = 0.15;% Shen uses 0.10*1e-2;
    E_max(j) = -log(f/(2*J+1))*pc.k*T/pc.icm_2_J;
    
    energies_g =(E_min(j):E_step:E_max(j)); % ground state energy

    %%%%%%%%%%%%%%%%%% find position R where potent curve height equal to energy value
    %     RR = zeros(length(energies_g),1);
    %     options = optimset('Diagnostics','off', 'Display','off');
    %     for i = 1: length(energies_g)
    %        fun = @(R) energies_g(i) - (U_X(R) + J*(J+1) * potJ(R));
    %   %          fun = @(R) energies_g(i) - (U_X(R) - J*(J+1)*potJ(R));
    %         %         RR(i) = vpasolve(energies_g(i) == U_X + J*(J+1)*potJ, RR);
    %         RR(i) = fsolve(fun,3.5,options); % RR       - where potent curve height equal to energy value
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %add R_far in the calculation
    [wvfns_g] = Free_Wvfns_NCV4(R,R_far,potpts_g_far(:,j),M_red,pc,energies_g); % R far for calculation wavefunction at the very far right end
    
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %     %Frank-Condon Factor, i stands For each ground state energy state,
% in each loop, calculate FCF between the ith ground state and numbound
% numbers of bound upper states

    for i = 1:length(energies_g)
        
        FCs_Mu = abs(wvfns_u.'*( murpts.*wvfns_g(:,i))).^2; % FCs_Mu each row corresponding to each upper vibrational state v'
        %       FCs_Mu = (abs((wvfns_u')*abs((wvfns_g(:,i)).*murpts))).^2; %WVNMs.^1;             % FC factors in a column vector
        S1 = FCs_Mu *(2*J+1)* exp(-energies_g(i)*pc.icm_2_J/(pc.k*T)); % Transition profile % coloumn vector
        E = energies_u-energies_g(i);           % wavelengths in wavenumber %column vector
        S2 = S1.*E; % Relative Intensity v*S(v), absorption coefficient k
        %         S2 = S1;
        ll = numbound;
        SS0((ll*(i-1)+1):ll*i,j) = S2;
        EE0((ll*(i-1)+1):ll*i,j) = E;
        
    end
    
    j=j+1;
end
    save('SS0','SS0');
    save('EE0','EE0');
if method ==1
    save('E_u','E_u');
    save('W_u','W_u');
elseif method ==2
    save('E_u_MN','E_u_MN'); %MN short for Matrix Numerov
    save('W_u_MN','W_u_MN');
end

SSS= reshape(SS0,[],1);
EEE= reshape(EE0,[],1);
% SSS(SSS==0)=[];
% EEE(EEE==0)=[];

end