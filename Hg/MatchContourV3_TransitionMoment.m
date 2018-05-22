%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gassian expansion has a better result then lorentizan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
% close all;

addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/Simulation')

M = csvread('HgSpectra300C10SPhotonControl.csv',0,0);
x = 1e7./M(:,1); %unit: nm
y = normc(M(:,2))*55;
figure;
plot (x,y,'g'); hold on;


%unit:Angstrom
% transition wavelength lambda to the ground state
lambda11 = 4865;
lambda12 = 3692;
lambda13 = 3749;

lambda(1,1:3)= [lambda11 lambda12 lambda13];

%unit:Angstrom
% transition wavelength lambda to the ground state
lambda21 = 3920;

lambda(1,4)= lambda21;

%unit:Angstrom
% transition wavelength lambda to the ground state
lambda31 = 3996;
lambda32 = 3897;
lambda33 = 3850;
lambda34 = 3936;
lambda35 = 3748;
lambda36 = 3738;
lambda37 = 2728;
lambda(1,5:11)= [lambda31 lambda32 lambda33 lambda34 lambda35 lambda36 lambda37];


p0 =1e8./lambda;

w = ones(1,11)*3500; %FWHM
%     w(1) = 2500;
%     w(4) = 4000;
%     w(8) = 3500;
  
    p = (1:1:1e3)*1e2;
    xx=zeros(11,1e3);
    yy=zeros(11,1e3);
    
k=1; %choose between Gassian or Lorentzian. k= 1 for Gaussian, else for Lorentzian
if k == 1;
    %--------------Gassian expansion---------------------------------

    for i = 1:11;
        xx(i,:) = (p0(i)-p)/(w(i)/2);
        yy(i,:) = 1/sqrt(2*pi)/(w(i)/2)*exp(-(p0(i)-p).^2/2/(w(i)/2)^2)*5e3;%Gaussian
        
    end
    
  
else
    %------------------Lorentizan--------------

    for i = 1:11;
        xx(i,:) = (p0(i)-p)/(w(i)/2);
        
        yy(i,:)= w(i)/2/pi./((w(i)/2)^2+(p0(i)-p).^2)*5e3; %Lorentian function
    end

end

%%%-----transition moment 
%---Dooh symmetrical linear configuration
yy(1,:) = yy(1,:)*.0550;%4865 angstrom
yy(2,:) = yy(2,:)*.108;%3692
yy(3,:) = yy(3,:)*.257;%3749
%---D3h equilateral triangular configuration
yy(4,:) = yy(4,:)*.102; %3920
%---C2v triangular configuration configuration
yy(5,:) = yy(5,:)*.0946;%3996
yy(6,:) = yy(6,:)*.0142;%3897
yy(7,:) = yy(7,:)*.0396;%3850
yy(8,:) = yy(8,:)*.110;%3936
yy(9,:) = yy(9,:)*.269;%3748
yy(10,:) = yy(10,:)*.04430;%3738    
yy(11,:) = yy(11,:)*1.343;%2728 

% plot(p,yy(1,:),'k.');hold on;
% plot(p,yy(9,:),'k.');hold on;
% plot(p,yy(11,:),'k.');hold on;

yy2 = yy(2:10,:);
yy3 = sum(yy2,1);
%%
yy4 = yy(1,:)*17.6 + yy3*3.08 + yy(11,:)*0.13; %add plots together

% yy5 = sum(yy,1);
plot(p,yy4,'r-');
xlabel('Wavenumber (cm^{-1})');
ylabel('Intensity (A.U.)');
xlim([1.2 4.2]*1e4);
ylim([0 4]);
% set(gca,'PlotBoxAspectRatio',[1.6 1 1]); %'DataAspectRatio',[1.6 1 1],...
title('under 248 nm excitation')
