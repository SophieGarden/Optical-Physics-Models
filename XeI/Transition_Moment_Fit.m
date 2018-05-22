clear;
clc;
% load('XeIShenExpData.mat');
% figure;
% Shen = sortrows([XeIShenExpData(:,1),XeIShenExpData(:,2)],1);
% ShenX=Shen(:,1);
% DIFShen=diff(ShenX);
% % plot(DIFShen);
% ShenY=Shen(:,2)./max(Shen(:,2));
% 
% [ShenX2,ia,ic] = unique(ShenX);
% DIFShen2=diff(ShenX2);
% % plot(DIFShen2);
% 
% 
% ShenY2=ShenY(ia);
% plot(ShenX2,ShenY2,'k'); %normalized exp data
% hold on;
% 
% [pks1,locs1] = findpeaks(ShenY2,ShenX2,'MinPeakDistance',80);
% findpeaks(ShenY2,ShenX2,'MinPeakDistance',80);
% 
% 
% load('ES2_ConstantTransitionMoment.mat')
% EE11= ES2_Mu1(:,1);
% n0=nnz(EE11);%# of zeros 
% n=length(EE11);
% EE22= ES2_Mu1(n-n0+1:end,1);
% SS22= ES2_Mu1(n-n0+1:end,2);
% plot(EE22, SS22*1.0,'-g');
% 
% [pks2,locs2] = findpeaks(SS22,EE22,'MinPeakDistance',80);
% findpeaks(SS22,EE22,'MinPeakDistance',80);


%%
% Pexp=locs1(14:14+33);
% Psim=locs2(11:11+33);
load('PeakExp');
load('PeakSim2');
Pexp=PeakExp(1,:);
Psim=PeakSim2(1,:);
%corresponding R values
% Rfit=@(x) 14.77564-6.06644E-4*x+8.30292E-9*x.^2;
Rfit=@(x) 6.34319-2.69317E-4*x+4.94456E-9*x.^2;
Rexp=Rfit(Pexp);
Rsim=Rfit(Psim);
% for ii=1:33
% V=EE22(Pexp(ii));
% end
Hexp=PeaksExp(:,2);
Hsim=PeaksSim(:,2);
MuV=(Hexp./Hsim).^(.5);
% MuV=MuV-min(MuV);
MuV=MuV./max(MuV);
R=2.2:5e-3:5.7;
% mur = @(R) 1.25^2./((1.25^2+(R-3.43).^2).*(1+exp(20*(R-4.3))))+.82^2./((.82^2+(R-3.60).^2).*(1+exp(-20*(R-4.30))));
% mur = @(R) 1^2./((1^2+(R-3.734).^2).*(1+exp(20*(R-4.00))))+.3479^2./((.3479^2+(R-3.8920).^2).*(1+exp(-20*(R-4.0))));
% mur = @(R) .4545^2./((.4545^2+(R-3.703).^2).*(1+exp(20*(R-4.00))))+.4187^2./((.4187^2+(R-3.771).^2).*(1+exp(-20*(R-4.0))));

%        a1 =     .5  ;% (0.4957, 0.5408)
%        b1 =       3.7;%  (3.691, 3.705)
%        a2 =      0.6 ;%  (0.5324, 0.7666)
%        b2 =       3.64 ;%  (3.478, 3.679)
% mur = @(R) a1^2./((a1^2+(R-b1).^2).*(1+exp(20*(R-4.4))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-20*(R-4.4))));

% five parameters
a1 =      0.4679 ;% (0.4493, 0.4866)
b1 =       3.707 ;% (3.699, 3.715)
d =       24.99  ;%(0.2788, 49.69)
c =       4.176  ;%(4.085, 4.268)
a2 =       0.208 ;% (0.1288, 0.2873)
b2 =       4.003 ;% (3.903, 4.102)
mur = @(R)  a1^2./((a1^2+(R-b1).^2).*(1+exp(d*(R-c))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-d*(R-c))));
murpts = mur(R);

 figure;
subplot(1,2,1)
plot(Pexp,MuV,'ro'); hold on;
xlabel('(Transition Energy (cm^{-1})');
ylabel('Normalized Transition Moment');

subplot(1,2,2)
plot(Rsim,MuV,'bo'); hold on; 
plot(R,murpts,'k-');hold on; 
xlabel('R (Angstrom)');
ylabel('Normalized Transition Moment');
title({'Dots are Calculated from peak height difference';' of Experiment and Simulation spectra.';'  Black curve is fit for blue dots, Red curve is adjusted fit'},'FontSize',14)

%% fit for transition moment

% % % %
% foexp = fitoptions('Method','NonlinearLeastSquares',...
%                 'Lower',[0 1 0 1 ],...
%                 'Upper',[1.25 10 10 10],...
%                 'StartPoint',[.4545 3.703 .4187 3.771]);
% Muexp = fittype('a1^2./((a1^2+(R-b1).^2).*(1+exp(20*(R-4.0))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-20*(R-4.0))))',...
%             'dependent',{'y'},'independent',{'R'},...
%             'coefficients',{'a1','b1','a2','b2'},'options',foexp);
%%%%
foexp = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0 1 .1 0 0 1 ],...
    'Upper',[2 10 100 10 2 10],...
    'StartPoint',[.4679 3.707 24.98 4.177 .208 4.033]);
Muexp = fittype('a1^2./((a1^2+(R-b1).^2).*(1+exp(d*(R-c))))+a2^2./((a2^2+(R-b2).^2).*(1+exp(-d*(R-c))))',...
    'dependent',{'y'},'independent',{'R'},...
    'coefficients',{'a1','b1','d','c','a2','b2'},'options',foexp);
% Addterm = fit(Rexp(3:end-5),MuV(3:end-5),Muexp)
Rexp2=zeros(length(Rexp),1);
for i=1:length(Rexp)
    if Rexp(i)<3.7
        Rexp2(i,1)=Rexp(i)-0.05;
    else
        Rexp2(i,1)=Rexp(i)+0.05;
    end
end

Addterm = fit(Rexp2,MuV,Muexp) %adjusted fit:fit for a widely distributed data
Addterm2 = fit(Rexp,MuV,Muexp) % best fit
Co = coeffvalues(Addterm);
plot(Addterm);
%%
% % % %4.066e4-4.467e4
% Pexp=locs1(14:14+33);
% Psim=locs2(11:11+33);
% Pdif=Pexp-Psim;
% % % %4.3063e4-4.467e4
% Pexp2=locs1(36:53);
% Psim2=locs2(31:48);
% Pdif2=Pexp2-Psim2;
%
% figure;
% plot(Psim,Pdif,'k*');hold on;
%
% % % %
% foexp = fitoptions('Method','NonlinearLeastSquares',...
%                 'Lower',[0 1e-5 0],...
%                 'Upper',[1 1 30],...
%                 'StartPoint',[2.527e-12 0.001 1]);
% UBexp = fittype('a*exp(x*k)+c',...
%             'dependent',{'y'},'independent',{'x'},...
%             'coefficients',{'a','k','c'},'options',foexp);
% % % %
% %         fopoly = fitoptions('Method','NonlinearLeastSquares',...
% %                 'Lower',[0 1e-5 0],...
% %                 'Upper',[1 1 30],...
% %                 'StartPoint',[2.527e-12 0.001 1]);
% % UBpoly = fittype('a*exp(x*k)+c',...
% %             'dependent',{'y'},'independent',{'x'},...
% %             'coefficients',{'a','k','c'},'options',fopoly);
%
% % Addterm = fit(Psim,Pdif,'poly9')
% %    Co = coeffvalues(Addterm);
% % plot(Addterm);
% figure;
% plot(Psim2,Pdif2,'k*');hold on;
% Addterm2 = fit(Psim2,Pdif2,'poly1')
%    Co2 = coeffvalues(Addterm2);
% plot(Addterm2);