%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 06/14/2016
%%%%
%%%% (1)Calulate and plot all the simulated waveforms for 335 nm and 485 nm
%%%%  band emission at One Temperatures set with TT
%%%% (2) compare calculated waveforms with experimentally recorded
%%%% waveforms of 335 nm and 485 nm.
%%%%
%%%%
%%%%
%%%%
%%%%
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
tic
global T n
%%----------------------------plot exp data at 266C------------------------
%addpath('HgCellwithVaccum/266nmPumped/SamePumpEnergy/266C');
order = 16;
if order ==1
    TT = 107;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\107C')
elseif order ==2
    TT = 138;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\138C')
elseif order ==3
    TT = 152;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\152C')
elseif order ==4
    TT = 162;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\162C')
elseif order ==5
    TT = 173;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\173C')
elseif order ==6
    TT = 182;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\182C')
elseif order ==7
    TT = 195;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\195C')
elseif order ==8
    TT = 204;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\204C')
elseif order ==9
    TT = 218;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\218C')
elseif order ==10
    TT = 229;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\229C')
elseif order ==11
    TT = 237;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\237C')
elseif order ==12
    TT = 266;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\266C')
elseif order ==13
    TT = 284;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\284C')
elseif order ==14
    TT = 298;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\298C')
elseif order ==15
    TT = 309;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\309C')
elseif order ==16
    TT = 319;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\319C')
elseif order ==17
    TT = 327;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\327C')
elseif order ==18
    TT = 335;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\335C')
elseif order ==19
    TT = 341;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\341C')
elseif order ==20
    TT = 348;
    addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\348C')
end




B335 = 7; %bandwidth of filter for 335 nm is 7 nm while for 485nm is 20 nm
B485 = 20; %bandwidth
magnify = 10; % 335 nm pump energy = 200uj; for 485 nm is 200uj for same energy case

M335 = csvread('335nm-200uj-2.csv',23,0);
x335 = M335(:,1)-1.3043e-8; %unit: second %correct the difference between photodiode and pmt
y335_pump = M335(:,3)-mean(M335(10:1400,3));%remove background noise
y335_pump2 = y335_pump./max(y335_pump);
y335_0 = - M335(:,2)-mean(-M335(10:1400,2));%remove background noise
Ratio = 1; %change to be same height with simulated curve
y335 =Ratio* y335_0./max(y335_pump);%normalization and magnification
y335 = y335/(B335/(B335+B485))/magnify;
%---filter
windowSize = 100;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y335_pump2 = filter(b,a,y335_pump2);
ratiooo=1;
y335 = filter(b,a,y335)*ratiooo;
% plot (x335,y335_pump2,'y.',x335,y335,'k-'); hold on;

M485 = csvread('485nm-20uj-2.csv',23,0);
x485 = M485(:,1)-1.3043e-8; %unit:  second
y485_pump = M485(:,3)-mean(M485(10:1400,3));%remove background noise
y485_pump2 = y485_pump./max(y485_pump);
y485_0 = - M485(:,2)-mean(-M485(10:1400,2));%remove background noise
y485 = Ratio* y485_0./max(y485_pump);%normalization and magnification
y485 = y485/(B485/(B335+B485));
%---filter
windowSize = 100;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y485_pump2 = filter(b,a,y485_pump2);
y485 = filter(b,a,y485);
%plot (x485,y485_pump2,'y.',x485,y485,'k-'); hold on;

%%-------------------------plotsimulation data at T=order---------------------
NA = 6.02214129e23;% mol-1, Avogadro constant
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

i=1;

Y2 = zeros(2e6,5);
T2 = zeros(2e6,1);
N1 = 0;

 T = TT + 273.15;
%-----------------
tau = 1-T/Tc;
P = Pc*exp((Tc/T)*(a1*tau+a2*tau^1.89+a3*tau^2+a4*tau^8+a5*tau^8.5+a6*tau^9));%Mpa
rho= P*1e6/(R*T);%m-3
n =  rho*NA*1e-6;%cm-3
N1(i) =n;
%-------------------
[T1,Y] = ode15s(@RateFunctionWendyV9FiveEqn,[0 2e-4]*1e6,[1 0 0 0 0]);
S = size(Y);
T2(1:S(1),i) = T1;
Y2(1:S(1),1:S(2),i) = Y;
% Y31(1:S(1),i) = Y(:,1);
% Y32(1:S(1),i) = Y(:,2);
% Y33(1:S(1),i) = Y(:,3);
% Y34(1:S(1),i) = Y(:,4);
% plot(T1,Y(:,1),'k.',T1,Y(:,2),'r.',T1,Y(:,3),'g.',T1,Y(:,3),'c.'); hold on;


T2(T2 == 0) = NaN;
Y2(Y2 == 0) = NaN;
T2 = T2./1e6;
N2 = transpose(N1);
csvwrite('HgsimulationN2_singleT.dat',N2);
%type HgsimulationN2.dat
csvwrite('HgsimulationY2_singleT.dat',Y2);
%type HgsimulationY2.dat
csvwrite('HgsimulationT2_singleT.dat',T2);
%type HgsimulationT2.dat
%%%

%%
j = 1;
Y3(:,1,j) = Y2(:,1,j);
Y3(:,2,j) = Y2(:,2,j); %*A335/(A335+A485);
Y4(:,2,j) =  Y3(:,2,j)/max( Y3(:,2,j))* max(y335(x335>10e-6)); %try to compare 

Y3(:,4,j) = Y2(:,4,j);%*A485/(A335+A485);
Y4(:,4,j) =  Y3(:,4,j)/max( Y3(:,4,j))*max(y485(x485>0.5e-6))*1;
%%
figure;  hold on;
subplot(2,1,1);
hold on;
plot (x335*1e6,y335_pump2,'c.',x335*1e6,y335,'k-');
plot(T2(:,j)*1e6,Y3(:,1,j),'g-');  %pump
plot(T2(:,j)*1e6,Y4(:,2,j),'b-.'); %335
axis([-1e-7*1e6 20e-5*1e6 -0.0002 0.02]);
xlabel('Time (micro second)');
ylabel('State Population');
title('335 nm');
legend('335Pump','335nm exp','N1 Population','335nm Simulation');

subplot(2,1,2);
hold on;
plot (x485*1e6,y485_pump2,'c.',x485*1e6,y485,'k-');
plot(T2(:,j)*1e6,Y3(:,1,j),'g-');  %pump
plot(T2(:,j)*1e6,Y2(:,2,j),'m:');  %state 2
plot(T2(:,j)*1e6,Y2(:,3,j),'m-');  %state 3n
plot(T2(:,j)*1e6,Y4(:,4,j),'r-.');  %485
axis([-1e-6*1e6  2e-4*1e6 -0.0002 0.2]);
xlabel('Time (micro second)');
ylabel('State Population');
title('485 nm');
legend('485Pump','485nm exp','N1 Population','N2 (335nm upper) Population','N3 Population','485nm Simulation');

% Y21=Y2(:,1,:);
%waterfalsave('pqfile.csv','T2','Y2','N2');l(N1,T2,Y31);
toc

%%
figure;  hold on;
subplot(2,1,1);
hold on;
plot (x335*1e6,y335,'k-');
% plot(T2(:,j)*1e6,Y3(:,1,j),'g-');  %pump
plot(T2(:,j)*1e6,Y4(:,2,j)*3,'b-.'); %335
axis([-1e-7*1e6 10e-5*1e6 -0.0002 0.02]);
xlabel('Time (micro second)');
ylabel('State Population');
title('335 nm');
legend('335 exp','335nm Simulation');

subplot(2,1,2);
hold on;
plot (x485*1e6,y485,'k-');
% plot(T2(:,j)*1e6,Y3(:,1,j),'g-');  %pump
% plot(T2(:,j)*1e6,Y2(:,2,j),'m:');  %state 2
% plot(T2(:,j)*1e6,Y2(:,3,j),'m-');  %state 3n
plot(T2(:,j)*1e6,Y4(:,4,j)*0.8,'r-.');  %485
axis([-1e-6*1e6  1e-4*1e6 -0.0002 0.2]);
xlabel('Time (micro second)');
ylabel('State Population');
title('485 nm');
legend('485 exp','485nm Simulation');

A1=horzcat(x335*1e6,y335,x485*1e6,y485);
csvwrite('319C_exp.dat',A1)
A2 = horzcat(T2(:,j)*1e6,Y4(:,2,j),Y4(:,4,j));
csvwrite('319C_simulation.dat',A2)
