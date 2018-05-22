clear;
close all;
i=1;
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

Y2 = zeros(2e5,4);
T2 = zeros(2e5,1);
N1 = zeros(1,20);
E23 = 2800; %cm-1 energy separation of dimer upper and middle level
E34 = 3800; %cm-1 energy separation of dimer middle and trimer upper level
E15 = 1690; %cm-1 energy separation of 3p1 and 3P0 level
E12 = 6900;
for T = [107,138,152,162,173,182,195,204,218,229,237,266,284,298,309,319,327,335,341,348,400] + 273.15
tau = 1-T/Tc;
    P = Pc*exp((Tc/T)*(a1*tau+a2*tau^1.89+a3*tau^2+a4*tau^8+a5*tau^8.5+a6*tau^9));%Mpa
    rho = P*1e6/(R*T);%m-3
    n =  rho*NA*1e-6;%cm-3
    N1(i) = n;
    k = 0.69503; %Blotzamann constant in unit of cm-1 k-1
kT = k*T; %cm-1
%plot(T,n,'r^-');hold on;
    %%
 rrr=1;%((204+273.15)/T)^1.4;
%rrr=1.7243e+04*(1/T)^1.4;
rrrr=rrr;
% plot(T,rrr,'r^-');hold on;
%%
k24 = 1e-34 * rrrr; %cm^6 sec-1, (between 1e-34 and 3e-31),trimer three body formation coefficient from dimer middle level 
r24 = k24 * n^2;

k12 = 5.35e-35; %Wendy's %fix
%r0 = k0 * n^2; %cm^6 sec^-1, dimer three body formation coefficient from 3P1;
r12 = k12 * n^2; %dimer three body formation rate from near 3P1;
r13 = r12/3*7; %/3 * 7; %dimer three body formation rate from near 3P1;

plot(n,r12,'r^-');hold on;
        xlabel('Hg Number Density (10^{18} cm^{-3})');
        ylabel('335 nm Fluoresence Decay Rate (10^{4} sec^{-1})');
k32 = 11.93e-14;
r32 = k32 * n;
plot(n,r32,'mo-');hold on;

k34 = 1.4e-34; %Wendy's %fix
r34 = k34 * n^2;
% plot(n,r34,'g^-');hold on;
% 
k43 = 2e-12 * 1; %0.001*rrr c; %cm^3 sec-1
r43 = k43 * exp(-E34/kT) * n;
% plot(n,r43,'b*-');hold on;

k23 = 10e-14; % Wendy's data between 1.75~2.03 e=15
r23 = k23 * n; 
% plot(n,r23,'b^-');hold on;

% k0 = 1.6e-31 * rrr; %1.6/1.8e-31;
% 
% r0 = k0 * n^2;
% r00 = 1e6;
% ratio= r00/r0;
% plot(T,ratio,'g^-');hold on;
%     T0 = 348+273.15;
%     T1 = 107+273.15;
%     cc3 = exp(-(T-T0)^2/ 3.0395e+06);

% T20 = 348 + 273.15;
% T1 = 107 + 273.15;
% ccc = log(30)/log(T1/T20);
% aaa = 1/T20^ccc;
% cc3 = aaa*T^ccc;
% Tc = T-273.15;
% %cc3 = exp(-(T-T0))^2/4.7800e+208;
%     
%       Y3 = cc3;
%       plot(Tc,Y3,'r^-');hold on;
%       
%     cc = 20;
%     Y = 4.2870 * (T/300)^(-2);
%     plot(Tc,Y,'bo-');hold on;
%    % R =8.3144621;
%    %Y2 = exp(+1450/R/T);
%    %plot(T,Y2,'k*-');hold on;
%      
%      Y1 = exp((204 + 273.15 - T)/20);
%      plot(Tc,Y1,'k*-');hold on;
%      %Y4 = (204 + 273.15 - T + 1) * 0.5;
%      %plot(Tc,Y4,'c*-');hold on;
%      
%      
%    		
% %138C-348C	
% 	y00=	-0.1165;
% A1=	1.2409E12;
% t1	=17.51415;
% %107C-348C	
% 	%Value	Standard Error
% % y00=1.29702;	%0.96081
% % A1=1.66003E16;	%8.88498E15
% % t1=12.48772;	%0.21957
% % a	=1.25418E16;%	6.18265E15
% % b	=0.92373;%	0.0012
% % rrr(i) = a*b^T	;
% 
% 
% %----------
% 
% A0=	-65.93163	;%183.23183
% A1=	0.8282;%	1.49835
% A2=	-0.00295;%	0.00456
% A3=	4.13882E-6;%;%	6.10248E-6
% A4=	-2.07362E-9;%	3.04075E-9
% 
% cccc = A0 + A1*T + A2*T^2 + A3*T^3 + A4*T^4;		
% rrr(i)=exp(cccc);
% %-----------------------------
% plot(Tc,rrr(i),'m*-');hold on;
% i=i+1;
% rrrt=transpose(rrr);
% 
% 
% y4 = A1*exp(-T/t1) + y00;	
% plot(Tc,y4,'c*-');hold on;
end