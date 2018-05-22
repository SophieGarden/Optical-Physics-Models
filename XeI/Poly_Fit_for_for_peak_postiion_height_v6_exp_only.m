
%% Find the peaks info of experimental data
load('XeIShenExpData.mat');
figure;
Shen = sortrows([XeIShenExpData(:,1),XeIShenExpData(:,2)],1);
ShenX = Shen(:,1);
% DIFShen=diff(ShenX);
% plot(DIFShen);
ShenY=Shen(:,2)./max(Shen(:,2));
[ShenX2,ia,~] = unique(ShenX);% there are replicate data in Shen
DIFShen2=diff(ShenX2);
% plot(DIFShen2);
ShenY2=ShenY(ia);

ShenY22=ShenY2;%save the orginal data in a different name vector named ShenY22 for later
%plot use

%--- remove noise spike
ShenY2(ShenX2>3.96065e+04 & ShenX2<3.9623e+04)=ShenY2(ShenX2>3.96065e+04 & ShenX2<3.9623e+04)-0.0163;
ShenY2(ShenX2>4.005912727e+04 & ShenX2<4.00764659e+04)=ShenY2(ShenX2>4.005912727e+04 & ShenX2<4.00764659e+04)- 0.0687;
ShenY2(ShenX2>3.985731083000000e+04 & ShenX2<3.989181340000000e+04)=ShenY2(ShenX2>3.985731083000000e+04 & ShenX2<3.989181340000000e+04)- 0.0082;

% ShenY2(249:259)=ShenY20(249:259)-0.02;
plot(ShenX2,ShenY2,'k'); hold on; %normalized exp data

%---remove baseline----------
ShenY2= msbackadj(ShenX2, ShenY2);
ShenY2= msbackadj(ShenX2, ShenY2);
% ShenY2= msbackadj(ShenX2, ShenY2);

plot(ShenX2,ShenY2,'r'); hold on; 

ExpData=[ShenX2, ShenY22,ShenY2];
% save('ExpData','ExpData');

[pks1,locs1] = findpeaks(ShenY2,ShenX2,'MinPeakDistance',80);
findpeaks(ShenY2,ShenX2,'MinPeakDistance',80);

% figure;
plot(ShenX2,-ShenY2,'r'); hold on; %normalized exp data
hold on;
[pks11,locs11] = findpeaks(-ShenY2,ShenX2,'MinPeakDistance',80);
findpeaks(-ShenY2,ShenX2,'MinPeakDistance',80);


N=49;
PeaksExpFit=zeros(N,2);
PeaksExpFit_trough=zeros(N,2);

load('PeakExp1.mat');
for i=1:N %1:51-73
    a=PeakExp1(3,i);
    if i<N
        b=PeakExp1(3,i+1);
    else
        b=PeakExp1(3,N)+80;
    end
    ShenX3=ShenX2(ShenX2>a & ShenX2<b);
    ShenY3=ShenY2(ShenX2>a & ShenX2<b);

    [p,~,mu] = polyfit(ShenX3,ShenY3,4);
    x1 = linspace(min(ShenX3),max(ShenX3));
    f1 = polyval(p,x1,[],mu);
    [P_fit,I_fit] = max(f1);
    
    %[P_fit,I_fit]  =findpeaks(f1);
    plot(x1(I_fit),P_fit,'r^'); hold on;
    PeaksExpFit(i,:)=[P_fit x1(I_fit)];
    plot(x1,f1,'r-','Linewidth',1.5);hold on;
    plot([x1(I_fit) x1(I_fit)],[-1 1],'c-');hold on;
end

for i=1:N %1:51-73 % get the bottom peak value
   if i==1
    a=PeakExp1(3,i);
   else
    a=PeaksExpFit(i-1,2);
   end
   b=PeaksExpFit(i,2);
    ShenX4=ShenX2(ShenX2>a & ShenX2<b);
    ShenY4=-ShenY2(ShenX2>a & ShenX2<b);
    %plot(ShenX4, ShenY4,'r-');hold on;
    [p2,~,mu2] = polyfit(ShenX4,ShenY4,4);
    x2 = linspace(min(ShenX4),max(ShenX4));
    f2 = polyval(p2,x2,[],mu2);
    [P_fit2,I_fit2]=max(f2);
    PeaksExpFit_trough(i,:)=[P_fit2 x2(I_fit2)];
    plot(x2,-f2,'b-');hold on;
    plot([x2(I_fit2) x2(I_fit2)],[-1 1],'m-');hold on;
end
title('Exoperimental');
%% %---experimental results:peak position and peak heights
PeakExp1(1,1:49)= PeaksExpFit(1:49,2);%peak position
PeakExp1(2,1:49) = PeaksExpFit(1:49,1)+PeaksExpFit_trough(1:49,1);%peak height
PeakExp1(3,1:49) = PeaksExpFit_trough(1:49,2);%trough position
PeakExp1(4,1:49) = PeaksExpFit(1:49,1);%peak max position
save('PeakExp1','PeakExp1');
% figure;
plot(PeakExp1(1,1:N),PeakExp1(4,1:N),'bo'); hold on;
