function dy = RateFunctionWendyV9FiveEqn(t,y)
%rate equation is constructed
global T n 
dy = zeros(5,1);    % a column vector
k = 0.69503; %Blotzamann constant in unit of cm-1 k-1
kT = k*T; %cm-1

E23 = 2800; %cm-1 energy separation of dimer upper and middle level
E34 = 3800; %cm-1 energy separation of dimer middle and trimer upper level
E15 = 1690; %cm-1 energy separation of 3p1 and 3P0 level
% E12 = 6900;
%----------------=============================================================
% a=[1.00000000000000e+16,0.0787451125322215,0.820224937072272];
% rrr= a(1)*exp(-a(2)*T)+a(3);
% radjust = 1;
% rrr=((266+273.15)/T)^1.4*radjust;
% rrr = 1;
% rrrr = rrr;


k12 = 1.65e-31; %stock CPL 67 1977
r12 = k12 * n^2; %dimer three body formation rate from near 3P1;

k53 = 1.55e-31; %stock CPL 67 1977
r53 = k53*n^2; %dimer three body formation rate from near 3P0;

k15 = 2.8e-13;%stock CPL 67 1977
r15 = k15 * n;

r51 = r15 * 3 *exp(-E15/kT);%stock CPL 67 1977

k24 = 5.4e-35; %cm^6 sec-1, Wendy,trimer three body formation coefficient from dimer middle level 
r24 = k24 * n^2;
%----------------
% A335 = 1e6; %spontaneous decay rate of lu, 335nm emission upper state
% A485 = 5.8e4; %A485~=A335/20,spontaneous decay rate of 485nm upper state

%set k23 according to k32
% k32 = 9.4e-14; %Wendy's
% %k23 = 1.6e-10; % * cm^3 sec^-1, collisional coefficient of dimer upper and middle state
% k23 = k32/exp(-E23/kT); % Wendy's data

%set k32 according to k23
k23 = 1.2e-15; % * cm^3 sec^-1, collisional coefficient of dimer upper and middle state
k32 = k23 *exp(-E23/kT); % Wendy's data

r23 = k23 * n;
r32 = k32 * n;

k34 = 1.8e-33*1; % Figen LOPE
r34 = k34 * n^2;

k4trimer=1.4e-34;%wendy's
r4trimer = k4trimer*n^2;

k423 = 5.9e-16; %k42+k43 %cm^3 sec-1; make sure that at largest[Hg], r43<r34
r423 = k423*n;

r42 = r423*0;
r43 = r423*1;

%spontaneous emission
A335 = 5.2e4; %wendy's A335 
A485 = 4.5e4;%Wendy's A485

%quenching loss/spontaneous emission
rq1 = 1.2e4; %stock CPL 67 1977 %1.2-2.5;
rq3 = 6.3e4;%Figen
rq5 = 1e4;%stock CPL 67 1977

rq2 = 0;%0.6e4;% add to fit better
rq4 = 1.2e4; % add to fit better
%diffusion loss
rd1 = 0;%10^20/n; %stock CPL 67 1977
rd5 = 0;%1.04e21/n; %stock CPL 67 1977
%%
P = 0; 
%P = rectangularPulse(0,1e-9,t); %this will make the calculation time longer
dy(1) = -(r12 + r15 + rq1 + rd1) * y(1) + r51 * y(5) + P; %3P1 atomic state population
dy(2) = r12 * y(1) + r32 * y(3) + r42 * y(4) - (A335 + r23 + r24 +rq2) * y(2); %dimer upper state
dy(3) = r53 * y(1) + r23 * y(2) + r43 * y(4) - (r32 + r34 + rq3) * y(3); %dimer middle state
dy(4) = r24 * y(2) + r34 * y(3) - (A485 + r423 + r4trimer +rq4) * y(4); %trimer upper state
dy(5) = r15 * y(1) - ( r51 + r53 + rq5 +rd5) * y(5);

dy = 1e-6*dy; %change time unit to [1/1e2 s]
