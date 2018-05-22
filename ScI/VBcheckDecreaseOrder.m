clear all;
clc;

%In order to get decreasing, we need We'<We'' ( WeXe'>WeXe'' )
%=D1 pia = B1 pia

%X1 sigma+, ground state
% Vel=0;
% Wel = 277.187; %ground state constants
% WeXel = 0.833;

%D1 pia
Veu = 15630.060;
Weu = 237.965;
WeXeu = 0.877;

%B1pia
Vel = 4409.708;
Wel = 243.432;
WeXel = 0.628;

% some pia
% Vel = 22285.31;
% Wel = 243.432;
% WeXel = 0.628;

deltaV1 = -1; %Vu-Vl=deltaV
deltaV2 = deltaV1+1;
m=0;%order of vibrational state of the UPPER electronic state of 1st group
n=0;%order of vibrational state of the UPPER electronic state of 2nd group
peakno = 10;

El = zeros(1,30);
Eu = zeros(1,30);
DE1 = zeros(1,peakno);
DE2 = zeros(1,peakno);
Ve = Veu - Vel;
for i=0:1:130
    v = i;
    El(i+1) = Wel*(v+0.5)-WeXel*(v+0.5)^2;
    Eu(i+1) = Ve + Weu*(v+0.5)-WeXeu*(v+0.5)^2;
end

%   for s=0:1:peakno+max(m,n)+k+15
%                 v = s;
%                 % ground state calculated energy value % El(1,2,3,...)
%                 % corresponding to V=0,1,2,...
%                 El(s+1,ii) = Welf(ii)*(v+0.5)-WeXelf(ii)*(v+0.5)^2;%Vel+
%                 Eu(s+1,ii) = Vef(ii)+Weuf(ii)*(v+0.5)-WeXeuf(ii)*(v+0.5)^2;%Vel+
%             end



for ss = 0:1:peakno-1
    
    DE1(ss+1) = -El(ss+m-deltaV1+1) + Eu(ss+1+m);
    DE2(ss+1) = -El(ss+n-deltaV2+1) + Eu(ss+1+n);
    
    %spacing between two groups
%     DE12(ss+1) = DE2(ss+1)-DE1(ss+1);
    
end

p(1) = -WeXeu+WeXel;
p(2) = Weu-Wel-WeXeu+WeXel*(-2*deltaV1+1);
p(3) = Ve+Weu/2-Wel*(-deltaV1+0.5)+WeXel*((-deltaV1+0.5)^2)-WeXeu/4;

%  p2(2) = -WeXeu2+WeXel;
%  p2(3) = Weu2-Wel-WeXeu2*(2*deltaV2+1)+WeXel;
%  p2(4) = Ve2+Weu2*(deltaV2+0.5)-WeXeu2*((deltaV2+0.5)^2)-Wel/2+WeXel/4;