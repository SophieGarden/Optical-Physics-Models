%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % generate spectra, the main code
% % % V8 find peaks, two pieces, first 10 peaks use large spacing
% % % V7 add continue to break out of the loop
% % % v6 simulation data,smooth for 1e4 data points, find trough first, then fit for peaks.
% % % v4 merge peaks diff code together to the main code
% % % v3, add iteration of the parameters, grid method. --6/6/2017
% % % add all the 16 parameters as variables.--05/28/2017 Wendy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
% close all;
% tic;
%experiment data plot
load('PeakExp1.mat'); %exp peak information % 4 rows
load('PeakSim1.mat'); %sim peak information
load('ExpData.mat'); %experimental curve
X1= ExpData(:,1); % wavenumber
Y1= ExpData(:,2); % original data
Y22= ExpData(:,3); % noise removed
Y2= msbackadj(X1, Y22); %baseline removed
Y2= msbackadj(X1, Y2);%baseline removed

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
%----plot index---
indicator_plot=1; % 1 for plot, 0 for no plot

% %peak positions:
N_first=1;
N_last=10;
N=48+1;
PeaksSimFit=zeros(N,2);
PeaksSimFit_trough=zeros(N,2);


                %begin calculations             
                %without transition moment
               [EE,SS] = Simulated_Spectra_WendyV13_XeI_function_No_TransitionMoment   (b_B,beta_B,De_l,beta_l,Re_l,De_r,beta_r,Re_r,a_off_l,Rs_off_l,a_on_l,Rs_on_l,a_off_r,Rs_off_r,a_on_r,Rs_on_r);
               %with transition moment   
%                [EE,SS] = Simulated_Spectra_WendyV13_XeI_function(b_B,beta_B,De_l,beta_l,Re_l,De_r,beta_r,Re_r,a_off_l,Rs_off_l,a_on_l,Rs_on_l,a_off_r,Rs_off_r,a_on_r,Rs_on_r);

                ES = sortrows([EE,SS],1);
                % plot(ES(:,1), ES(:,2),'-b'); hold on;
                Wvl = 1e7./ES(:,1); %wavelength in [nm]
                % adding the effect of refractie index has no effect
                n = 1 ; %+ 0.05792105/(238.0185-Wvl.^(-2)) + 0.00167917/(57.362-Wvl.^(-2)); % air refractive index
                SS1 = ES(:,2);
                SS1 = SS1/max(SS1);
                ES1 = [Wvl,SS1];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Parameters for Jstep=1
                WvlRel = 2e-3*0.5*0.8; %nm
                instrRes = 0.06*50*0.5;  %nm  %assumed the resolution of the monochrometer used in Gary's 1993 paper since it is not provided
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % %Parameters for Jstep=5
%                 WvlRel = 2e-3*0.5*4; %nm
%                 instrRes = 0.06*50*3; %nm  %assumed the resolution of the monochrometer used in Gary's 1993 paper since it is not provided
                %  %%%%%%%%%%%%%%%%%%%%%%%%%
                %Parameters for Jstep=10
                % WvlRel = 2e-3*0.5*4; %nm
                % instrRes = 0.06*50*2; %nm  %assumed the resolution of the monochrometer used in Gary's 1993 paper since it is not provided
                
                %Parameters for Jstep=20
                % WvlRel = 2e-3*0.5*8; %nm
                % instrRes = 0.06*50*2;  %nm  %assumed the resolution of the monochrometer used in Gary's 1993 paper since it is not provided
                %
                %Parameters for SingleJ
                % WvlRel = 2e-3*0.5*20; %nm
                % instrRes = 0.06*50;  %nm  %assumed the resolution of the monochrometer used in Gary's 1993 paper since it is not provided
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %P, Q, R brach
                % The actual FWHM of spectral lines is the convolution of the intrinsic line-width by this instrumental line-width.
                instr_wvl = -5*instrRes:WvlRel:5*instrRes;
                filter = instrRes/2/pi./((instrRes/2)^2+instr_wvl.^2); %Lorentian function
                % filter = 1/instrRes/sqrt(pi).*exp(-instr_wvl.^2./(instrRes/2)^2); %Gaussian function
                SS2 = conv(ES1(:,2),filter,'same');
                %%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!
                SS2 = SS2/max(SS2);
                EE2 = 1e7./ES1(:,1).*n';
                ES_sim=[EE2, SS2];
                EE11= ES_sim(:,1);
                n0=nnz(EE11); %# of zeros
                n=length(EE11);
                EE22= ES_sim(n-n0+1:end,1);
                SS22= ES_sim(n-n0+1:end,2);
                [EE22,ia,ic] = unique(EE22);% remove replicate data
                SS22=SS22(ia);
                % % % %                 ES2 = [EE22, SS22];
                % % % %                 save('ES2','ES2');
                %% % % % % plot sim and exp spectra
                if indicator_plot==1
                    figure;
                    plot(X1,Y1,'k'); %normalized exp data
                    hold on;
                    plot(X1,Y2,'b','linewidth',1.0);
                    hold on;
                    plot(EE22, SS22,'-c');
                    hold on;
                    
                    legend('Experimental data','Simulation')
                    xlim([3.91e4 4.52e4]);%for XeI
                    ylim([-.05 1.05]);
                    xlabel('Dye Laser Frequency (cm^{-1})');
                    ylabel ('B-X and B-A Transition Relative Intensity');
                    title('XeI Photoassociative Excitation Spectra');
                    set(gca,'FontSize',20)
                end
                %% compare sim and exp difference, peak positions and peak heights
                %% %%%%%(step 1)-----------find troughs of simulation data------------%%%%%%%%%%%%%%%
                % PeakSim1=PeakExp1;
                EE22=EE22(EE22>3.90e4 & EE22<4.53e4);
                SS22=SS22(EE22>3.90e4 & EE22<4.53e4);
                SS22=SS22./max(SS22);
                %---remove baseline---
                SS22b= msbackadj(EE22, SS22);%plot(EE22, SS22,'-b','LineWidth',2.0);hold on;
                SS22b= msbackadj(EE22, SS22b);%plot(EE22, SS22,'-r','LineWidth',1.5);hold on;
                SS222= smooth(SS22b,1e4,'moving');
                EE2p1= EE22(EE22>3.96e4&EE22<4.121e4);
                SS2p1= -SS222(EE22>3.96e4&EE22<4.121e4);
                EE2p2= EE22(EE22>4.121e4&EE22<4.51e4);
                SS2p2= -SS222(EE22>4.121e4&EE22<4.51e4);
                
             
                space1=110;
                [pks_t_1,locs_t_1] = findpeaks(SS2p1,EE2p1,'MinPeakDistance',space1);
                space2=70;
                [pks_t_2,locs_t_2] = findpeaks(SS2p2,EE2p2,'MinPeakDistance',space2);
                
                if indicator_plot==1
                    
                    plot(EE22, SS22b,'-c','LineWidth',2.0);hold on;
                    % find peaks only plot their own region, i.e. half of
                    % the whole graph for each command
%                     findpeaks(SS2p1,EE2p1,'MinPeakDistance',space1);hold on;%plot out
%                     findpeaks(SS2p2,EE2p2,'MinPeakDistance',space2);hold on;%plot out
                end
                
                 locs_t=[locs_t_1;locs_t_2];
                 pks_t=[pks_t_1;pks_t_2];
             
               
                %% %%%%%%%%(step 2)---------find peaks of simulation data
                
                for i=1:N
                    if i==1
                        a=locs_t(1)-200;
                    else
                        a=locs_t(i-1);
                    end
                    b=locs_t(i);
                    
                    EE23=EE22(EE22>=a & EE22<=b);
                    SS23=SS22b(EE22>=a & EE22<=b);
                    [p,~,mu] = polyfit(EE23,SS23,4);
                    x_local = linspace(min(EE23),max(EE23));
                    fun2 = polyval(p,x_local,[],mu);
                    [P_fit_sim,I_fit_sim]=max(fun2);
                    
                    PeaksSimFit(i,:)=[P_fit_sim x_local(I_fit_sim)];
                    if indicator_plot==1
                        plot(x_local,fun2,'r-','linewidth',1);hold on;
                        plot([x_local(I_fit_sim) x_local(I_fit_sim)],[-1 1],'r-','linewidth',1.5);hold on;
                    end
                    
                end
                
                for i=1:N %1:51-73 % get the bottom trough value
                    if i==1
                        a=PeakExp1(3,i);
                    else
                        a=PeaksSimFit(i-1,2);
                    end
                    b=PeaksSimFit(i,2);
                                        
                    EE24=EE22(EE22>a & EE22<b);
                    SS24=-SS22b(EE22>a & EE22<b); % reversed
                    %plot(ShenX4, ShenY4,'r-');hold on;
                    [p2,~,mu2] = polyfit(EE24,SS24,4);
                    x2_local = linspace(min(EE24),max(EE24));
                    fun2 = polyval(p2,x2_local,[],mu2);
                    [P_fit2,I_fit2]= max(fun2);
                    PeaksSimFit_trough(i,:)=[P_fit2 x2_local(I_fit2)];
                    if indicator_plot==1
                        plot(x2_local,-fun2,'b-');hold on;
                        plot([x2_local(I_fit2) x2_local(I_fit2)],[-1 1],'b-');hold on;
                    end
                end
                title('Simulation');
                % %--Simulated results:peak position and peak heights
                PeakSim2(1,1:N) = PeaksSimFit(1:N,2);%peak position
                PeakSim2(2,1:N) = PeaksSimFit(1:N,1)+PeaksSimFit_trough(1:N,1);%peak height
                PeakSim2(3,1:N) = PeaksSimFit_trough(1:N,2);%trough position
                PeakSim2(4,1:N) = PeaksSimFit(1:N,1);%peak max position
                % save('PeakSim2','PeakSim2');% does not need to be saved
                
                % --difference between experiment and simulation
                % --difference add together, similar to standard deviation
                
                %PeakExp1(1,1)=3.961161490959596e+04;
                diff_p=(PeakSim2(1,1:N)-PeakExp1(1,1:N)).^2/N; % position difference
                diff_position=sqrt(sum(diff_p,2));
                diff_h=(PeakSim2(2,1:N)./PeakExp1(2,1:N)-mean(PeakSim2(2,1:N)./PeakExp1(2,1:N))).^2/N; % height ratio VS 1 %%improve!!!!
                diff_height=sqrt(sum(diff_h,2));
                %% save the data
                   DataTDMConstant=[EE22, SS22];
                    
                    csvwrite('DataTDMConstant.dat', DataTDMConstant);
%%
% % % % figure;
% % % % subplot(2,1,1)
% % % % plot((1:1:N),PeakSim2(1,1:N)-PeakExp1(1,1:N),'bo');
% % % % xlabel('Peak Number');
% % % % ylabel('Position difference (Sim-Exp)');
% % % % subplot(2,1,2)
% % % % plot((1:1:N),PeakSim2(2,1:N)./PeakExp1(2,1:N),'ro');
% % % % xlabel('Peak Number');
% % % % ylabel('Height difference (Sim/Exp)');

%% 
figure;
%---sim---
h=SS22-SS22b;
plot(EE22,SS22,'r-','LineWidth',1.5); hold on;
plot(EE22,h,'k-',EE22,SS22b,'b-',EE22,SS22b+h,'c-'); hold on;
EE_h=EE22(EE22>3.94e4 & EE22<4.51e4);
hh=h(EE22>3.94e4 & EE22<4.51e4);
[p_h,~,mu_h] = polyfit(EE_h,hh,21);
x_local_h = linspace(min(EE_h),max(EE_h));
fun_h = polyval(p_h,x_local_h,[],mu_h);
plot(x_local_h,fun_h,'r-','linewidth',1);hold on;

delta_h=polyval(p_h,PeakSim2(1,1:N),[],mu_h);

Isim_y=delta_h+PeakSim2(4,1:N);
plot(PeakSim2(1,1:N),Isim_y,'ro');

%% --exp--
figure;
h2=Y1-Y2;
plot(X1,Y1,'r-','LineWidth',1.5); hold on;
plot(X1,h2,'k-',X1,Y2,'b-',X1,Y2+h2,'c-'); hold on;

for i=1:N
    %max difference between two points is 4.1502
X_h=X1(X1>=(PeakExp1(1,i) - 2.1) & X1<=(PeakExp1(1,i) + 2.1));
I_h1=find(X1>=(PeakExp1(1,i) - 2.1) & X1<=(PeakExp1(1,i) + 2.1));
I_h=round(mean(I_h1));
delta_h2(1,i)=h2(I_h);
end
Iexp_y=delta_h2+PeakExp1(4,1:N);

plot(PeakExp1(1,1:N),PeakExp1(4,1:N),'bo');
plot(PeakExp1(1,1:N),Iexp_y,'ro');


%% %corresponding R values
% Rfit=@(x) 14.77564-6.06644E-4*x+8.30292E-9*x.^2;
% Rfit=@(x) 6.34319-2.69317E-4*x+4.94456E-9*x.^2;
figure;
R=3.17:5e-4:5.7; %3.1660 is the monoton point
U_B = @(R) 73160.89 + b_B*exp(-beta_B*R)-1.16141e5*R.^(-1)-3.9e4*R.^(-3)-5.3e5*R.^(-4)-6.3e5*R.^(-6);
%---U_X-----
U_L_Morse = @(R) De_l*((exp(-beta_l*(R-Re_l))-1).^2-1);
U_R_Morse = @(R) De_r*((exp(-beta_r*(R-Re_r))-1).^2-1); 

U_X = @(R) U_L_Morse(R)./(1+exp(a_off_l*(R-Rs_off_l)))...
    +U_R_Morse(R)./((1+exp(-a_on_l*(R-Rs_on_l))).*(1+exp(a_off_r*(R-Rs_off_r))))...
-(1.73e6*R.^(-6)+1.21e7*R.^(-8))./(1+exp(-a_on_r*(R-Rs_on_r))); 

 UBR=U_B(R)';
UXR=U_X(R)';

Diff_BX=U_B(R)-U_X(R);
Diff_BX=Diff_BX(Diff_BX>3.91e4 & Diff_BX < 4.52e4);
R=R(Diff_BX>3.91e4 & Diff_BX < 4.52e4);
plot(Diff_BX, R,'k-');hold on;


[p_bx,~,mu_bx] =polyfit(Diff_BX,R,10);

x_local_bx = linspace(min(Diff_BX),max(Diff_BX));
fun_bx = polyval(p_bx,x_local_bx,[],mu_bx);
plot(x_local_bx,fun_bx,'r-','linewidth',1);hold on;

% load('PeakExpT.mat');
Pexp=PeakExp1(1,1:N);
% Rexp=Rfit(Pexp);
% Rsim=Rfit(PeakSim2(1,1:N));
Rexp=polyval(p_bx,Pexp,[],mu_bx);
plot(Pexp, Rexp,'b*');hold on;

%%
figure;
Hexp=Iexp_y;%PeakExpT(2,1:N);
Hsim=Isim_y;
Mue=(Hexp./Hsim).^(.5);
% try use baseline difference
%Mue=(delta_h2./delta_h).^(.5);
Mue=Mue./max(Mue);

%manually change points for test
%Mue=[0.552501266250982,0.781509705296526,0.872535245977445,0.907214592473497,0.906665597988582,0.927064052115872,0.942030213108040,0.960013930853821,0.970429309031394,0.984546380019311,0.988918269455341,0.992881861618215,1,0.981116400058566,0.991325819288339,0.971482126614088,0.978913158221505,0.986044300339691,0.970263521366561,0.963863402064888,0.961536021917643,0.961648494551448,0.945497391634188,0.939836129962768,0.917136341329774,0.898028309460871,0.866492253064717,0.846469138167896,0.808302300966393,0.779325573770731,0.763766072094645,0.747085828128481,0.735268234151824,0.717403430982489,0.706144706158549,0.693972401798910,0.671226726557594,0.655361362952682,0.642990485612099,0.616716080951868,0.609505640016902,0.592734515809802,0.577366107844568,0.580597291590633,0.550442245974027,0.519647490419341,0.489631250200784,0.478727153367988,0.441847440114590];

n=7;
p = polyfit(Rexp',Mue',n);
mur0 = @(x)  p(1)*x.^n+p(2)*x.^(n-1)+p(3)*x.^(n-2)+p(4)*x.^(n-3)+p(5)*x.^(n-4)+p(6)*x.^(n-5)+p(7)*x.^(n-6)+p(8)*x.^(n-7);
% %  
% n=6;
% p = polyfit(Rexp',Mue',n);
%  mur0 = @(x)  p(1)*x.^n+p(2)*x.^(n-1)+p(3)*x.^(n-2)+p(4)*x.^(n-3)+p(5)*x.^(n-4)+p(6)*x.^(n-5)+p(7)*x.^(n-6);
 
%  n=5;
% p = polyfit(Rexp',Mue',n);
%  mur0 = @(x)  p(1)*x.^n+p(2)*x.^(n-1)+p(3)*x.^(n-2)+p(4)*x.^(n-3)+p(5)*x.^(n-4)+p(6)*x.^(n-5);
%  
 
%  n=3;
% p = polyfit(Rexp',Mue',n);
%  mur0 = @(x)  p(1)*x.^n+p(2)*x.^(n-1)+p(3)*x.^(n-2)+p(4)*x.^(n-3);
 
murpts0 = mur0(R);

% fit for transition moment

%====================================
%---plot

% five parameters
   a1 =      0.6466 ;% (0.6046, 0.6887)
       b1 =       3.704 ;% (3.684, 3.725)
       d =       4.651  ;%(1.906, 7.396)
       c =       4.656  ;%(4.353, 4.96)
       a2 =     0.06297;%  (-0.09029, 0.2162)
       b2 =       3.888 ;% (3.822, 3.954)

mur = @(R)  a1^2./((a1^2+(R-b1).^2).*(1+exp(d*(R-c))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-d*(R-c))));
murpts = mur(R);

 figure;
subplot(1,2,1)
plot(Pexp,Mue,'ro'); hold on;
xlabel('(Transition Energy (cm^{-1})');
ylabel('Normalized Transition Moment');

subplot(1,2,2)
 plot(Rexp,Mue,'bo'); hold on; 
plot(R,murpts,'k-');hold on; 
plot(R,murpts0,'r-');hold on; 
xlabel('R (Angstrom)');
ylabel('Normalized Transition Moment');
title({'Dots are Calculated from peak height difference';' of Experiment and Simulation spectra.';'  Black curve is fit for blue dots, Red curve is adjusted fit'},'FontSize',14)
