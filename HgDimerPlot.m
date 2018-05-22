clear;
close all;
clc;
% mass= 200.59*1.66053904e-27;
% u = mass/2;%reduced mass
% beta20=sqrt(8*3.141592653^2*3e8*0.4250*u/(6.626070040e-34)); 
% %unit: (sqrt(m*s-1*cm-1*kg/(kg·m2/s2*s))=sqrt(cm^-1*m^-1)=cm^-1*0.1
% beta2=beta20*10;
% in the order of X D A
De = [379.5  8229.3 1.279874558896620e+04]; %cm^-1
Re =[3.605 2.591 2.6023];  %a^-1
beta=[1.2377  1.7921 1.54782071642272]; %a^-1
% Umorse=zeros(10,2);
% r =zeros(10,2);

for j=1:3
    i=1;
    for R= 2:0.01:7
       Umorse(i,j) = De(j).*(1-exp(-beta(j)*(R-Re(j)))).^2;
        r(i)=R;
        i=i+1;
    end
    
end
r=transpose(r);
Ux0 = Umorse(:,1);
UD0 = Umorse(:,2);%+39412.3-De(2)
Ux = Ux0/8065.73;
UD = UD0/8065.73 + 4.1114;
plot(r,Ux,'r-',r,UD,'b-');hold on;

UA0 = Umorse(:,3);%+28800
UA = UA0/8065.73 + 3.6393- 0.0917;% manullay add a difference 0.0917 to force the two states merge together
plot(r,UA,'g-');hold on;
xlabel('Internuclear distance (Angstrom)');
ylabel('Energy(ev)');
%=========================
mass= 200.59*1.66053904e-27;
u = mass/2;%reduced mass
h = 6.626070040e-34;
c=3e8;
reA = 2.6023;%angstrom
%DeA = 2.91e-9; %cm^-1
wexeA = 0.4024; %cm^-1
weA = 143.530; %cm-1
BeA = sqrt(8*pi^2*c*u*wexeA/h)*1e-9;%ansgtrom^-1
DeA= weA^2/(4*wexeA);

