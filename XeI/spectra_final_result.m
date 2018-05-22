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
tic;
%experiment data plot
load('PeakExp1.mat'); %exp peak information
load('PeakSim1.mat'); %sim peak information
load('ExpData.mat'); %experimental curve
X1= ExpData(:,1); % wavenumber
Y1= ExpData(:,2); % original data
Y2= ExpData(:,3); % baseline removed

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
N_first=3;
N_last=48;
N=48;
PeaksSimFit=zeros(N,2);
PeaksSimFit_trough=zeros(N,2);


                %begin calculations             
                %withOUT transition moment
     %           [EE,SS] = Simulated_Spectra_WendyV13_XeI_function_No_TransitionMoment   (b_B,beta_B,De_l,beta_l,Re_l,De_r,beta_r,Re_r,a_off_l,Rs_off_l,a_on_l,Rs_on_l,a_off_r,Rs_off_r,a_on_r,Rs_on_r);
               %with transition moment   
               [EE,SS] = Simulated_Spectra_WendyV13_XeI_function(b_B,beta_B,De_l,beta_l,Re_l,De_r,beta_r,Re_r,a_off_l,Rs_off_l,a_on_l,Rs_on_l,a_off_r,Rs_off_r,a_on_r,Rs_on_r);

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
                
                %nomalization in the plot region
                EE22=EE22(EE22>3.90e4 & EE22<4.55e4);
                SS22=SS22(EE22>3.90e4 & EE22<4.55e4);
                SS22=SS22./max(SS22);
                
                % % % %                 ES2 = [EE22, SS22];
                % % % %                 save('ES2','ES2');
                %% % % % % plot sim and exp spectra
                if indicator_plot==1
                    figure;
                    plot(X1,Y1,'k'); %normalized exp data
                    hold on;
                    plot(X1,Y2,'k','linewidth',1.0);
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
                
%                 %nomalization in the plot region
%                 EE22=EE22(EE22>3.90e4 & EE22<4.53e4);
%                 SS22=SS22(EE22>3.90e4 & EE22<4.53e4);
%                 SS22=SS22./max(SS22);
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
                % save('PeakSim2','PeakSim2');% does not need to be saved
                
                % --difference between experiment and simulation
                % --difference add together, similar to standard deviation
        %%     
                PeakExp1(1,1)=3.961161490959596e+04;
                diff_p=PeakSim2(1,1:N)-PeakExp1(1,1:N); % position difference
                diff_p_2=diff_p-mean(diff_p);
                diff_position=sqrt(sum(diff_p_2.^2,2)/N);%standard deviation
                 diff_h=PeakSim2(2,1:N)./PeakExp1(2,1:N); % height ratio VS 1 %%improve!!!!
                 diff_h_2=diff_h-mean(diff_h);
   %             diff_h=(PeakSim2(2,1:N)./PeakExp1(2,1:N)-1).^2; % height ratio VS 1 %%improve!!!!

                diff_height=sqrt(sum(diff_h_2.^2,2)/N);
                %% Part of the peaks, not all %%N_first, N_last
                %%%% %--Simulated results:peak position and peak heights
                PeakSim2(1,N_first:N_last) = PeaksSimFit(N_first:N_last,2);%peak position
                PeakSim2(2,N_first:N_last) = PeaksSimFit(N_first:N_last,1)+PeaksSimFit_trough(N_first:N_last,1);%peak height
                PeakSim2(3,N_first:N_last) = PeaksSimFit_trough(N_first:N_last,2);%trough position
                
                %% --difference between experiment and simulation
                % --difference add together, similar to standard deviation
                
                                PeakExp1(1,1)=3.961161490959596e+04;
                diff_p_part=PeakSim2(1,N_first:N_last)-PeakExp1(1,N_first:N_last); % position difference
                diff_p_2_part=diff_p_part-mean(diff_p_part);
                diff_position_part=sqrt(sum(diff_p_2_part.^2,2)/(N_last-N_first+1));%standard deviation
                 diff_h_part=PeakSim2(2,N_first:N_last)./PeakExp1(2,N_first:N_last); % height ratio VS 1 %%improve!!!!
                 diff_h_2_part=diff_h_part-mean(diff_h_part);
                  diff_height_part=sqrt(sum(diff_h_2_part.^2,2)/(N_last-N_first+1));
                
%                 diff_p_part=(PeakSim2(1,N_first:N_last)-PeakExp1(1,N_first:N_last)).^2/(N_last-N_first+1); % position difference
%                 diff_position_part=sqrt(sum(diff_p_part,2));
%                 diff_h_part=(PeakSim2(2,N_first:N_last)./PeakExp1(2,N_first:N_last)-mean(PeakSim2(2,N_first:N_last)./PeakExp1(2,N_first:N_last))).^2/(N_last-N_first+1); % height ratio VS 1 %%improve!!!!
%                 diff_height_part=sqrt(sum(diff_h_part,2));
                
%%
figure;
subplot(2,1,1)
plot((1:1:48),PeakSim2(1,1:N)-PeakExp1(1,1:N),'bo');
xlabel('Peak Number');
ylabel({'Position difference';' Sim-Exp (cm^{-1})'});
subplot(2,1,2)
plot((1:1:48),PeakSim2(2,1:N)./PeakExp1(2,1:N)-1,'ro');
xlabel('Peak Number');
ylabel({'Height difference Ratio';'Sim/Exp - 1'});
%%
figure;

                    plot(X1,Y1,'b','LineWidth',1.5); %normalized exp data
                    hold on;
                    


                    plot(EE22, SS22,'-r','LineWidth',1.5);
                    hold on;
                    DataSimulation=[EE22, SS22];
                    %save('DataSimulation','DataSimulation');
                    csvwrite('DataSimulation.dat', DataSimulation);
                    legend('Experimental data','Simulation')
                    xlim([3.91e4 4.5e4]);%for XeI
                    ylim([-.05 1.05]);
                    xlabel('Dye Laser Frequency (cm^{-1})');
                    ylabel ('B-X and B-A Transition Relative Intensity');
                    title('XeI Photoassociative Excitation Spectra');
                    set(gca,'FontSize',20)
               
toc;
disp(toc/60)