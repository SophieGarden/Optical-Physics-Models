%function changing_reflectivity
%*************written by Wenting Chen****************
clear all;
clc;%���µ�λ�����ر�˵����Ϊ���ʵ�λ��
tic
global h c  Ner Nyb sgm56p sgm65p sgm12 sgm21 sgm56 sgm65 ... 
   lams1 lams2 alpha1 alpha2 alphap Ts1 Ts2 Tp N2 N3 N6 N6np

h=6.628e-34;%plunk constant
c0=2.9979245e8;%light velocity in vaccum
k=1.3806505e-23;%Boltzmann constant
pia=3.1416;%Բ����

t21=10e-3;%Er���ϼ����ܼ�����
t32=1e-7;%Er��N3�ܼ�����
t65=1.3e-3;%Yb���ܼ�����

r1=30e-6/2;%core radius=Diameter/2
r2=400e-6/2;%inner cladding radius=Diameter-short axis of D shape inner cladding/2
r3=640e-6/2;%m. polymer,outer cladding radius
Aco=pia*r1^2;%core area
Acl=pia*r2^2;%innner cladding area
nco=1.449;%refract index of core
%ncl=1.432;%refract index of inner cladding
c=c0/nco;
NNNN=0;
for R1er=0.01:-0.005:0.005
    NNNN=NNNN+1;
%----------�߲���-----------
Ner=5.1e25;  %4e25; %%����Ũ��Er  ����ʵ�ʲ��������������û����ô�࣡����
Nyb=7.6e26; %5e26;  %%����Ũ��Yb (Nyb/Ner=14) the ratio of Nyb/Ner is approximately 4~20/8~30

C2=3.5e-24+(Ner-4.4e25)*2.4107e-49;%������ת��ϵ�� Ner>4.4e25�������
R61=1e-22+(Nyb-1e25)*3.96e-49;%energy transfer coefficient������  1e25<Nyb<1e27  
R35=R61;%2.371e-22;%��������ת��ϵ��2.371e-22   
%-----------�Ͳ���-------------------
%Ner=3.7e25;   %3.7e25;      %2.450e25;  %4e25; %%����Ũ��Er 
%Nyb=5.3e26;   %5.3e26;    %3.526e26; %5e26;  %%����Ũ��Yb (Nyb/Ner=14) the ratio of Nyb/Ner is approximately 4~20/8~30

%C2=3e-24; %1e-23;%3e-24;%������ת��ϵ��
%R61=2.371e-22; %3.5e-22;%2.371e-22;%energy transfer coefficient 
%R35=R61;%2.371e-22;%��������ת��ϵ��2.371e-22   
%------fnp-----------------------
fnp=0.15;  %fnp���Լ���Ϊ a simple function of NYb/ NEr
Nybp=Nyb*(1-fnp);%����Er Yb �������ݵ�
Nybnp=Nyb*fnp;%������Nybnp
%---------------------------------------------------------------------------------------------------------
%R1er=R1er;%reflecty of fiber side-Er 1.00
R2er=R1er;%0.001;%%reflecty of the other side-Er  0.036
R1yb=R1er;%reflecty of fiber side-Yb 1.00
R2yb=R1er;%0.001;%%reflecty of fiber the other side-Yb  0.034

alpha1=3.2e-3;   %20db/km =0.004605  %3.2e-3;%1540nm�����
alpha2=5e-3;  %70db/km =0.0161  %5e-3;%1064nm�����
alphap=6.5e-3;  %70db/km =0.0161  %6.5e-3;%975nm�����

%�±� s1=1.5��m s2=1��m P=pump 975nm
lams1=1.543e-6;%�źŹⲨ��1
lams2=1.058e-6;%�źŹⲨ��2
lamp=0.975e-6;%���ֲ���

NAcore=0.22; %0.22; %0.16;%sqrt(ncore^2-ninner^2); % numerical aperture of the core
Vs1=2*pia*r1*NAcore/lams1;% V is the normalized parameter of the core
Vs2=2*pia*r1*NAcore/lams2;
Vp=2*pia*r1*NAcore/lamp;
ws1=r1*(0.65+1.619/Vs1^1.5+2.879/Vs1^6);% mode radius 
ws2=r1*(0.65+1.619/Vs2^1.5+2.879/Vs2^6);
wp=r1*(0.65+1.619/Vp^1.5+2.879/Vp^6);
Ts1=1-exp(-2*r1^2/ws1^2);% wavelength dependent overlap factor between the light field and the doped core Լ����1
Ts2=1-exp(-2*r1^2/ws2^2);%Լ����1
Tp=(1-exp(-2*r1^2/wp^2))*Aco/Acl;%Լ����0.005 %��Ȼ�Ƚ�С�����Ǳ��ֹ⻹�Ǽ�������������
%=============================================================================
L=5; %length of fiber
m=50; %Max ���� ----L/m �����������С����������ʾ���̲���� %�������M()��������


%------����洢�ռ�--------%�����ʱ�� u0 ������P56�ĵ���һ�£�����������������
%-----------PP_g(6,m+1,u,g), Ti(g,u)...sgm21(g,u) 
g0=3;
u0=10;
ZZ_g=zeros(1,m+1,u0,g0);
pp_g=zeros(6,m+1,u0,g0);
Ti=zeros(g0,u0);
Pump=zeros(g0,u0);
Per=zeros(g0,u0);
Pyb=zeros(g0,u0);
    
Sgm56p=zeros(g0,u0);
Sgm65p=zeros(g0,u0);
Sgm56=zeros(g0,u0);
Sgm65=zeros(g0,u0);
Sgm13p=zeros(g0,u0);
Sgm13=zeros(g0,u0);
Sgm12=zeros(g0,u0);
Sgm21=zeros(g0,u0);

%------�ֲ���--------------------------------------------------------------
% j=1��ʾ���£�j=2��ʾg�¶��£�g��ʾ�����¶�
M=zeros(1,2); %�����¶�
Zf=zeros(6,2); %6���ܼ���2���¶�

Fer=zeros(8,2,3,3);% (8,j,3,g)
Fyb=zeros(4,2,2,3);%(4,j,2,g)
%============================================================================
%-----------------�����¶ȶ����շ���ϵ����Ӱ��---------------------------------
%���¸�E�ܼ���λΪcm-1
E=zeros(6,8); %6���ܼ�
E(1,:)=[0 47 111 125 164 288 306 318]; %Er��̬�����ܼ� 8��  -n1
%E(1,:)=[0 26 51 59 125 133 201 168]; %Er��̬�����ܼ� 8��  -n1
E(2,:)=[6529 6563 6613 6641 6679 6737 6743 0]; %Er��һ����̬�����ܼ� 7�� -n2
E(3,:)=[10166 10186 10220 10256 10267 0 0 0]; %Er�ĵڶ�����̬���ܼ� 5�� -n3
E(4,:)=[0 0 0 0 0 0 0 0]; % �ò��� -n4
E(5,:)=[0 492 970 1365 0 0 0 0]; %Yb�Ļ�̬���ܼ� 4�� -n5
E(6,:)=[10239 10909 11689 0 0 0 0 0]; %Yb�ĵ�һ����̬���ܼ� 3��  -n6

Tr=300;%Tr=room temperature/����
M(1)=-h*c*100/k/Tr;%exp()�����h*c/��/k/T��cm-1�����m-1,�ʳ���100 %����

g=0; %g�������ݶȣ�T0��Χ100��
for Tair=273:50:273  %k temperature of the air  %ǰ����㲨�������ֲ��õ�Tr��room temperature
    g=g+1;
    u=0;
 for P56=20:50:470  %Launched pump power in one end. totoal=one*2
   u=u+1;
%------����洢�ռ�-------------   
p0=zeros(6,m+1);
p=zeros(6,m+1);Z=zeros(1,m+1);
p01=zeros(6,m+1);
p1=zeros(6,m+1);Z1=zeros(1,m+1);
p02=zeros(6,m+1);
p2=zeros(6,m+1);Z2=zeros(1,m+1);
p03=zeros(6,m+1);
p3=zeros(6,m+1);Z3=zeros(1,m+1);
eta=zeros(1,u);
%------����quantum defect ����Ĺ���core�¶�����---------------------
if u==1  %u==1��ʱ��Ҫ��ybû����ʱ��ʼ����Ϊ��ʱ����eta=0.4,��yb���������
 eta(1,u)=0.4; %��,���ֹ�ת�����ȵ�Ч�ʣ�Ҳ����δ�����յı��ֹ����
else
    eta(1,u)=1-(Per(g,u-1)+Pyb(g,u-1))/Pump(g,u-1); %��ǰ��һ�����eta������һ���㣬�򲽳�ԽС��etaֵԽ��
end
Pdiss=(P56*2)*eta(1,u);% Pdiss is mainly produced by the quantum defect, which is the absorbed pump power minus the output signal power at both 1&1.5 ��m
Q=Pdiss/(L*Aco);% Q = Pdiss/(L*Aco) is the heat density generated in the core of the fibre
I0=(1-eta(1,u))*(P56*2)/(L*2*pia*r2); %the luminescence intensity at r=r2;
alphaS=10; %m-1 the absorption coefficient in the wavelength range of spontaneous emission in the polymer coating
x=alphaS*I0*r2*exp(alphaS*r2);

ks=1.38;%Wm-1k-1  the thermal conductivity of the silica/quartz glass %��boltzmann Сдk����
kp=0.24;%Wm-1k-1   the thermal conductivity of the polymer outer cladding
ka=0.026;%Wm-1k-1   the thermal conductivity of the air
%H=25;%h=115W m-2 K-1 heat transfer/convective coefficient of the polymer-air interface(���¶����еģ����¶������½�,��������ȡһ�̶�ֵ��%��дH��plunk constantСдh����
gg=9.8;%m*s-2
Pr=0.71;
dA=1.2;%kg m-3
miuA=1.85e-5;%kg m-1 s-1
epson=0.85;
sigmab=5.669e-8;
%----------------boundary condition: continuous& Newton's cooling law------
syms cc2 cc5 cc6 cc7 cc8 nu hh%����Constant CC/cc ��֮ǰ�Ĺ���c,c0�Լ�cp conversion C2����

Gr=gg*(dA/miuA)^2*(2*r3)^3/Tair*(cc7+cc8*log(r3)-x*real(expint(alphaS*r3))/(alphaS*kp)-Tair);
%Nu=2.0/lambertw(0, 25/(2*(Gr*Pr)^(1/3)));%solve ��������̸�eq6=nu*exp(-2/nu)-0.16*(Gr*Pr)^(1/3);
%Hc=0.5*nu*ka/r3;
%Hr=epson*sigmab*((cc7+cc8*log(r3)-x*real(expint(alphaS*r3))/(alphaS*kp))^4-(Tair)^4)/((cc7+cc8*log(r3)-x*real(expint(alphaS*r3))/(alphaS*kp))-(Tair));
hh=0.5*nu*ka/r3+epson*sigmab*((cc7+cc8*log(r3)-x*real(expint(alphaS*r3))/(alphaS*kp))^4-(Tair)^4)/((cc7+cc8*log(r3)-x*real(expint(alphaS*r3))/(alphaS*kp))-(Tair));
eq1=cc2-Q*r1^2/(4*ks)-cc5-cc6*log(r1);
eq2=cc5+cc6*log(r2)-cc7-cc8*log(r2)+x*real(expint(alphaS*r2))/(alphaS*kp);
eq3=-Q*r1/(2*ks)-cc6/r1;
eq4=ks*cc6/r2-kp*(cc8/r2+x*exp(-alphaS*r2)/(alphaS*kp*r2));
eq5=cc8/r3+x*exp(-alphaS*r3)/(alphaS*kp*r3)-hh/kp*(Tair-cc7-cc8*log(r3)+x*real(expint(alphaS*r3))/(alphaS*kp));
eq6=nu*exp(-2/nu)-0.16*(Gr*Pr)^(1/3);

[CC2 CC5 CC6 CC7 CC8 Nu]=solve(eq1,eq2,eq3,eq4,eq5,eq6,'cc2','cc5','cc6','cc7','cc8','nu');
CC2=double(CC2);CC5=double(CC5);CC6=double(CC6);CC7=double(CC7);CC8=double(CC8);Nu=double(Nu);
Hc=0.5*Nu*ka/r3;
Hr=epson*sigmab*((CC7+CC8*log(r3)-x*real(expint(alphaS*r3))/(alphaS*kp))^4-(Tair)^4)/((CC7+CC8*log(r3)-x*real(expint(alphaS*r3))/(alphaS*kp))-(Tair));

T=CC2;       % center temperature of the core %T0�ǹ������������¶�, %melting temprature for Sio2:1982K;
T123=zeros(1,101);
R=zeros(1,101);
if g==1&&P56==470
for nr=1:1:101
    R(nr)=(nr-1)*r3/100;
    if R(nr)<r1||R(nr)==r1
        T123(nr)=CC2-Q*R(nr)^2/(4*ks); %kת�������϶ȡ�
    elseif R(nr)>r1&&(R(nr)<r2||R(nr)==r2)
         T123(nr)=CC5+CC6*log(R(nr));
    elseif R(nr)>r2
         T123(nr)=CC7+CC8*log(R(nr))-x*real(expint(alphaS*R(nr)))/(alphaS*kp);
    end
end
subplot(2,4,4);
RR=R*1e6; %��λ��m����ɦ�m
TT123=T123-273; %��λ��k����������¶�
plot(RR,TT123,'m-');hold on;
title('Temperature profile in the radical direction');
xlabel('Radius/��m');
ylabel('Temperature/��');
end
   %T2=T0+Q*r1^2/2/ks*(log(r2/r1)+ks/r2/H+0.5);% center temperature of the core %T0�ǹ������������¶�
   %Tcl=T0+Q*a^2/(2*b*H); %cladding temperature
   %T=Tc;%T=Tc-Q*a^2/12/K;%average temperature of the core of fiber %average temperature T is approximately Tc ʵ��core�¶ȿ��ܻ�ܸ�

M(2)=-h*c*100/k/T;%average temperature of core  
for j=1:2 %1=300k,2=���Ⱥ��ҿ���fiber core�¶����ߺ�Ӱ���core�¶�
%---------Zf��j), partition function �ֱ��upper �� lower levels ���루�Ƕ�upper/lower��������)--%�����Ĺ��˳���Z��������----------------------------------------------------------------------
Zf(1,j)=(exp(M(j)*E(1,1))+exp(M(j)*E(1,2))+exp(M(j)*E(1,3))+exp(M(j)*E(1,4))+exp(M(j)*E(1,5))+exp(M(j)*E(1,6))+exp(M(j)*E(1,7))+exp(M(j)*E(1,8)))/exp(M(j)*E(1,1));
Zf(2,j)=(exp(M(j)*E(2,1))+exp(M(j)*E(2,2))+exp(M(j)*E(2,3))+exp(M(j)*E(2,4))+exp(M(j)*E(2,5))+exp(M(j)*E(2,6))+exp(M(j)*E(2,7)))/exp(M(j)*E(2,1));
Zf(3,j)=(exp(M(j)*E(3,1))+exp(M(j)*E(3,2))+exp(M(j)*E(3,3))+exp(M(j)*E(3,4))+exp(M(j)*E(3,5)))/exp(M(j)*E(3,1));
Zf(5,j)=(exp(M(j)*E(5,1))+exp(M(j)*E(5,2))+exp(M(j)*E(5,3))+exp(M(j)*E(5,4)))/exp(M(j)*E(5,1));
Zf(6,j)=(exp(M(j)*E(6,1))+exp(M(j)*E(6,2))+exp(M(j)*E(6,3)))/exp(M(j)*E(6,1));
%--------pupulation %------------------------------------
Fer(1,j,1,g)=exp(M(j)*E(1,1))/Zf(1,j)/exp(M(j)*E(1,1));
Fer(2,j,1,g)=exp(M(j)*E(1,2))/Zf(1,j)/exp(M(j)*E(1,1));
Fer(3,j,1,g)=exp(M(j)*E(1,3))/Zf(1,j)/exp(M(j)*E(1,1));
Fer(4,j,1,g)=exp(M(j)*E(1,4))/Zf(1,j)/exp(M(j)*E(1,1));
Fer(5,j,1,g)=exp(M(j)*E(1,5))/Zf(1,j)/exp(M(j)*E(1,1));
Fer(6,j,1,g)=exp(M(j)*E(1,6))/Zf(1,j)/exp(M(j)*E(1,1));
Fer(7,j,1,g)=exp(M(j)*E(1,7))/Zf(1,j)/exp(M(j)*E(1,1));
Fer(8,j,1,g)=exp(M(j)*E(1,8))/Zf(1,j)/exp(M(j)*E(1,1));

Fer(1,j,2,g)=exp(M(j)*E(2,1))/Zf(2,j)/exp(M(j)*E(2,1));
Fer(2,j,2,g)=exp(M(j)*E(2,2))/Zf(2,j)/exp(M(j)*E(2,1));
Fer(3,j,2,g)=exp(M(j)*E(2,3))/Zf(2,j)/exp(M(j)*E(2,1));
Fer(4,j,2,g)=exp(M(j)*E(2,4))/Zf(2,j)/exp(M(j)*E(2,1));
Fer(5,j,2,g)=exp(M(j)*E(2,5))/Zf(2,j)/exp(M(j)*E(2,1));
Fer(6,j,2,g)=exp(M(j)*E(2,6))/Zf(2,j)/exp(M(j)*E(2,1));
Fer(7,j,2,g)=exp(M(j)*E(2,7))/Zf(2,j)/exp(M(j)*E(2,1));

Fer(1,j,3,g)=exp(M(j)*E(3,1))/Zf(3,j)/exp(M(j)*E(3,1));
Fer(2,j,3,g)=exp(M(j)*E(3,2))/Zf(3,j)/exp(M(j)*E(3,1));
Fer(3,j,3,g)=exp(M(j)*E(3,2))/Zf(3,j)/exp(M(j)*E(3,1));
Fer(4,j,3,g)=exp(M(j)*E(3,4))/Zf(3,j)/exp(M(j)*E(3,1));
Fer(5,j,3,g)=exp(M(j)*E(3,5))/Zf(3,j)/exp(M(j)*E(3,1));

Fyb(1,j,1,g)=exp(M(j)*E(5,1))/Zf(5,j)/exp(M(j)*E(5,1));
Fyb(2,j,1,g)=exp(M(j)*E(5,2))/Zf(5,j)/exp(M(j)*E(5,1));
Fyb(3,j,1,g)=exp(M(j)*E(5,3))/Zf(5,j)/exp(M(j)*E(5,1));
Fyb(4,j,1,g)=exp(M(j)*E(5,4))/Zf(5,j)/exp(M(j)*E(5,1));

Fyb(1,j,2,g)=exp(M(j)*E(6,1))/Zf(6,j)/exp(M(j)*E(6,1));
Fyb(2,j,2,g)=exp(M(j)*E(6,2))/Zf(6,j)/exp(M(j)*E(6,1));
Fyb(3,j,2,g)=exp(M(j)*E(6,3))/Zf(6,j)/exp(M(j)*E(6,1));
end
if g==3&&P56==470
subplot(2,4,1)
bar(Fer(:,:,1,g),'grouped');
colormap summer
%set(gca,'XTickLabel',{'a','b','c','d','e','f','g','h'})

legend('T ��','T ��');
title 'Er Stark sublevels �ֲ��������¶����߱仯����'
xlabel('Er Stark sublevels of 4I15/2');
ylabel('Occupation Ratio');
        
subplot(2,4,2)
bar(Fer(:,:,2,g),'grouped');
%set(gca,'XTickLabel',{'i','j','k','l','m','n','o'})

legend('T ��','T ��');
title 'Er Stark sublevels �ֲ��������¶����߱仯����'
xlabel('Er Stark sublevels E2');
ylabel('Occupation Ratio');

subplot(2,4,3)
bar(Fer(:,:,3,g),'grouped')
%set(gca,'XTickLabel',{'p','q','r','s','t'})
legend('T ��','T ��');
title 'Er Stark sublevels �ֲ��������¶����߱仯����'
xlabel('Er Stark sublevels E3');
ylabel('Occupation Ratio');
        
subplot(2,4,5)
bar(Fyb(:,:,1,g),'grouped')
%set(gca,'XTickLabel',{'a','b','c','d'})

legend('T ��','T ��');
title 'Yb Stark sublevels �ֲ��������¶����߱仯����'
xlabel('Yb Stark sublevels E5');
ylabel('Occupation Ratio');

subplot(2,4,6)
bar(Fyb(:,:,2,g),'grouped')
%set(gca,'XTickLabel',{'e','f','g'})

legend('T ��','T ��');
title 'Yb Stark sublevels �ֲ��������¶����߱仯����'
xlabel('Yb Stark sublevels E6');
ylabel('Occupation Ratio');
end
%---------------r: room temperature �µķ������--------------
sgm56pr=14e-25; %2.85e-24;  %14e-25;  %pump 975
sgm65pr=7e-25; %2.62e-24;  %7e-25;

sgm13pr=2e-25;
sgm13r=0;%�ο����ײ���������ĵ�transition cross section (m^2)���ǳ���(��������Ϊ��300k)�µĲ���ֵ��������ԣ���ĸ�ǣ��ϡ����ܼ����������� �������sublevel/���ܼ��ġ� e.g. ���ܼ���absorption cross section=absorption coefficient/(fl*N)

sgm56r=0;%1064
sgm65r=0.1e-25;%1064

sgm21r=5e-25; %5e-25;%2.5e-25;%1540
sgm12r=5e-25; %5e-25;%2e-25;
%---------------temperature of core �µ����շ������-------------------
Sgm65p(g,u)=sgm65pr*Zf(6,1)/Zf(6,2);%�����¶�Ӱ���� % E(6,1)->E(5,1);
Sgm56p(g,u)=Sgm65p(g,u)*Zf(6,2)/Zf(5,2)*exp(0); %Ezl=hv
%Sgm56p(g,u)=sgm56pr*Zf(5,1)/Zf(5,2)*exp((E(6,1)-E(5,1)-(E(6,1)-E(5,1)))*(M(2)-M(1)));

Sgm13p(g,u)=sgm13pr*Zf(1,1)/Zf(1,2)*exp((E(3,1)-E(1,1)-(E(3,4)-E(1,1)))*(M(2)-M(1)));%(E(3,4)<-E(1,1)
Sgm13(g,u)=0;%sgm13r=0;

Sgm65(g,u)=sgm65r*Zf(6,1)/Zf(6,2); %E(6,1)->E(5,2)
Sgm56(g,u)=0; %Attention��������Ϊsgm56<<sgm65,����Բ���һ��Ҫ��sgm56=0;��������Ts2*(n6*sgm65-n5*sgm56)-alpha2 < 0���������p(3,:)�����󷴶���С�������� %Sgm56(g,u)=Sgm65(g,u)*Zf(6,2)/Zf(5,2)*exp((E(6,1)-E(5,1)-(E(6,1)-E(5,2)))*M(2)); 

Sgm21(g,u)=sgm21r*Zf(2,1)/Zf(2,2);%E(2,2)->E(1,3)  ~1.55��m
Sgm12(g,u)=Sgm21(g,u)*Zf(2,2)/Zf(1,2)*exp((E(2,1)-E(1,1)-(E(2,2)-E(1,3)))*M(2));
%Sgm12(g,u)=sgm12r*Zf(1,1)/Zf(1,2)*exp((E(2,1)-E(1,1)-(E(2,2)-E(1,3)))*(M(2)-M(1))); 

%-----------------------
sgm65p=Sgm65p(g,u);
sgm56p=Sgm56p(g,u);
sgm13p=Sgm13p(g,u);
sgm13=Sgm13(g,u);
sgm65=Sgm65(g,u);
sgm56=Sgm56(g,u);
sgm21=Sgm21(g,u);
sgm12=Sgm12(g,u);
%==========================================================================

%----------------shooting method-------------------------------------------
delta1=100;
delta2=100;
delta3=100;
%--------give/guess the initial values,p0(5,1)/p0(2,1),p0(4,1),p0(6,1)--------
% v=1; %z=0�����������

p0(5,1)=P56;           %pump+0
p6L=P56;               %pump-L

p0(2,1)=p0(5,1)*0.4;            %ps1-0
p0(4,1)=p0(5,1)*0.2;             %ps2-0
p0(6,1)=1e-7*P56; %5e-5;             %pump-0

p0(1,1)=R1er*p0(2,1);      %ps1+0
p0(3,1)=R1yb*p0(4,1);      %ps2+0
%---Ԥ���ֵ���--------------------------------------------------------------------------------------------------------------------
%p0(:,1)=[10.3639555482135;6.72484593551156;1.35520348325411;12.1428289188094;16.2477923861786;7.44716556712002e-05];

while (abs(delta1)>2||abs(delta2)>2||abs(delta3)>2)      %����ֵ ������delta1,delta2>2, delta3>0.1

%---------------calculate the other 2 initial values,p10,p30----------------------------------
%v=1; %z=0�����������


%-----------------p0(v)��L�ĵ�v�εĳ�ֵ,v=1��Z=0~1��һ������˵ĳ�ֵ��p0(2)�ǵڶ��εĳ�ֵ-----------------
%-----------------p(1)��v=1��Z=0~1��һ������˵ĳ�ֵ��p(2)����L��Z=0~1��һ���Ҷ˵ļ���ֵ��Ҳ�ǵڶ��εĳ�ֵ-----p(v+1)�ǵ�v�ε��Ҷ�ֵ----------------
Z(1)=0; 
p(:,1)=p0(:,1);
%----------------ode45��M������L���ȵĹ��˷ָ���M����һ��һ������һ��Z=Zv+1���ĳ�ֵ--------------

for v=2:m+1
%---------------��Z=z0v���� ���ܼ��������ֲ� n1 n2 n3 n5 n6------------------    
 W12=Ts1*sgm12*(p0(1,v-1)+p0(2,v-1))*lams1/(h*c*Aco);
 W21=Ts1*sgm21*(p0(1,v-1)+p0(2,v-1))*lams1/(h*c*Aco);
 W13=(Tp*sgm13p*(p0(5,v-1)+p0(6,v-1))*lamp+Ts2*sgm13*(p0(3,v-1)+p0(4,v-1))*lams2)/(h*c*Aco);
 W56=(Tp*sgm56p*(p0(5,v-1)+p0(6,v-1))*lamp+Ts2*sgm56*(p0(3,v-1)+p0(4,v-1))*lams2)/(h*c*Aco);
 W65=(Tp*sgm65p*(p0(5,v-1)+p0(6,v-1))*lamp+Ts2*sgm65*(p0(3,v-1)+p0(4,v-1))*lams2)/(h*c*Aco);
 
syms n2 n3 n6
eq1=-n2/t21+n3/t32+W12*(Ner-n2-n3)-W21*n2-2*C2*n2^2;%dN2/dt=0
eq2=-n3/t32+W13*(Ner-n2-n3)+R61*n6*(Ner-n2-n3)-R35*n3*(Nybp-n6)+C2*n2^2;%dN3/dt=0
eq3=-n6/t65+W56*(Nybp-n6)-W65*n6-R61*n6*(Ner-n2-n3)+R35*n3*(Nybp-n6);%dN6/dt=0

[N2 N3 N6]=solve(eq1,eq2,eq3,'n2','n3','n6');
N2=double(N2);
N3=double(N3);
N6=double(N6);
for i=1:3
    if N2(i)<Ner&&N2(i)>0&&N6(i)>0 %�������⣬ֻ��һ������0��������
        N2=N2(i);
        N3=N3(i);
        N6=N6(i);
        break
    end
end
%---nonparticipatory Yb3+ -ions--------------------------------------------
syms n6np
eq4=-n6np/t65+W56*(Nybnp-n6np)-W65*n6np;%dN6p/dt=0

[N6np]=solve(eq4,'n6np');
N6np=double(N6np);

%----------------ode45 ��� ���䷽��----------------------------------------
z0v=[L/m*(v-2) L/m*(v-1)];
p0v=p0(:,v-1);
[z,y]=ode45(@divide_double_pumped_fnp,z0v,p0v);

x=length(z(:,1));
Z(v)=z(x); 
p(:,v)=y(x,:);

p0(:,v)=p(:,v);
end
%**************************************************************************
%##########################################################################
%--------------------------F1 F2 F3 ��Z=L���ı߽�����-----------------------
%F1 F2 F3�������Ų�Ҫ��Ŷ����ΪҪ���������б�ʵ�ʱ��ͳһ��
F1=p(2,m+1)-R2er*p(1,m+1);
F2=p(4,m+1)-R2yb*p(3,m+1);
F3=p(6,m+1)-p6L;
%##########################################################################

%-------------------------У�� ������ֵp0(2,1)��p0(4,1),p0(6,1)-------------
%----1.55��m & 1��m & 975nm ��Ӧ�ĳ�ֵ p20 p40 p60����б�ʵ�С��-------------
eps1=p0(2,1)/20;  %p0(2,1)/5;
eps2=p0(4,1)/20;   %p0(4,1)/5; %0.1;
eps3=p0(6,1)/20;   %p0(6,1)/5; %1e-6;

%------------------------p0(2,1)+eps1-----begin----------------------------
%**************************************************************************
%----------------guess the initial values----------------------------------
% v=1; %z=0�����������
p01(5,1)=p0(5,1);           %pump+0
p01(6,1)=p0(6,1);           %pump-0
p01(2,1)=p0(2,1)+eps1;        %ps1-0
p01(4,1)=p0(4,1);             %ps2-0
p01(1,1)=R1er*p01(2,1);      %ps1+0
p01(3,1)=R1yb*p01(4,1);      %ps2+0
%-----------------p0(v)��L�ĵ�v�εĳ�ֵ,v=1��Z=0~1��һ������˵ĳ�ֵ��p0(2)�ǵڶ��εĳ�ֵ-----------------
%-----------------p(1)��v=1��Z=0~1��һ������˵ĳ�ֵ��p(2)����L��Z=0~1��һ���Ҷ˵ļ���ֵ��Ҳ�ǵڶ��εĳ�ֵ-----p(v+1)�ǵ�v�ε��Ҷ�ֵ----------------
Z1(1)=0; 
p1(:,1)=p01(:,1);
%----------------ode45��M������L���ȵĹ��˷ָ���M����һ��һ������һ��Z=Zv+1���ĳ�ֵ--------------

for v=2:m+1
%---------------��Z=z0v���� ���ܼ��������ֲ� n1 n2 n3 n5 n6------------------    
 W12=Ts1*sgm12*(p01(1,v-1)+p01(2,v-1))*lams1/(h*c*Aco);
 W21=Ts1*sgm21*(p01(1,v-1)+p01(2,v-1))*lams1/(h*c*Aco);
 W13=(Tp*sgm13p*(p01(5,v-1)+p01(6,v-1))*lamp+Ts2*sgm13*(p01(3,v-1)+p01(4,v-1))*lams2)/(h*c*Aco);
 W56=(Tp*sgm56p*(p01(5,v-1)+p01(6,v-1))*lamp+Ts2*sgm56*(p01(3,v-1)+p01(4,v-1))*lams2)/(h*c*Aco);
 W65=(Tp*sgm65p*(p01(5,v-1)+p01(6,v-1))*lamp+Ts2*sgm65*(p01(3,v-1)+p01(4,v-1))*lams2)/(h*c*Aco);
 
syms n2 n3 n6
eq1=-n2/t21+n3/t32+W12*(Ner-n2-n3)-W21*n2-2*C2*n2^2;%dN2/dt=0
eq2=-n3/t32+W13*(Ner-n2-n3)+R61*n6*(Ner-n2-n3)-R35*n3*(Nybp-n6)+C2*n2^2;%dN3/dt=0
eq3=-n6/t65+W56*(Nybp-n6)-W65*n6-R61*n6*(Ner-n2-n3)+R35*n3*(Nybp-n6);%dN6/dt=0

[N2 N3 N6]=solve(eq1,eq2,eq3,'n2','n3','n6');
N2=double(N2);
N3=double(N3);
N6=double(N6);
for i=1:3
    if N2(i)<Ner&&N2(i)>0&&N6(i)>0 %�������⣬ֻ��һ������0��������
        N2=N2(i);
        N3=N3(i);
        N6=N6(i);
        break
    end
end
%---nonparticipatory Yb3+ -ions--------------------------------------------
syms n6np
eq4=-n6np/t65+W56*(Nybnp-n6np)-W65*n6np;%dN6p/dt=0

[N6np]=solve(eq4,'n6np');
N6np=double(N6np);
%----------------ode45 ��� ���䷽��----------------------------------------
z01v=[L/m*(v-2) L/m*(v-1)];
p01v=p01(:,v-1);
[z1,y1]=ode45(@divide_double_pumped_fnp,z01v,p01v);

x1=length(z1(:,1));
Z1(v)=z1(x1); 
p1(:,v)=y1(x1,:);

p01(:,v)=p1(:,v);

end
%**************************************************************************
%------------------------p0(2,1)+eps1-----end------------------------------

%------------------------p0(4,1)+eps2-----begin----------------------------
%**************************************************************************
%----------------guess the initial values----------------------------------
% v=1; %z=0�����������
p02(5,1)=p0(5,1);           %pump+0
p02(6,1)=p0(6,1);           %pump-0
p02(2,1)=p0(2,1);        %ps1-0
p02(4,1)=p0(4,1)+eps2;             %ps2-0
p02(1,1)=R1er*p02(2,1);      %ps1+0
p02(3,1)=R1yb*p02(4,1);      %ps2+0
%-----------------p0(v)��L�ĵ�v�εĳ�ֵ,v=1��Z=0~1��һ������˵ĳ�ֵ��p0(2)�ǵڶ��εĳ�ֵ-----------------
%-----------------p(1)��v=1��Z=0~1��һ������˵ĳ�ֵ��p(2)����L��Z=0~1��һ���Ҷ˵ļ���ֵ��Ҳ�ǵڶ��εĳ�ֵ-----p(v+1)�ǵ�v�ε��Ҷ�ֵ----------------
Z2(1)=0; 
p2(:,1)=p02(:,1);
%----------------ode45��M������L���ȵĹ��˷ָ���M����һ��һ������һ��Z=Zv+1���ĳ�ֵ--------------

for v=2:m+1
%---------------��Z=z0v���� ���ܼ��������ֲ� n1 n2 n3 n5 n6------------------    
 W12=Ts1*sgm12*(p02(1,v-1)+p02(2,v-1))*lams1/(h*c*Aco);
 W21=Ts1*sgm21*(p02(1,v-1)+p02(2,v-1))*lams1/(h*c*Aco);
 W13=(Tp*sgm13p*(p02(5,v-1)+p02(6,v-1))*lamp+Ts2*sgm13*(p02(3,v-1)+p02(4,v-1))*lams2)/(h*c*Aco);
 W56=(Tp*sgm56p*(p02(5,v-1)+p02(6,v-1))*lamp+Ts2*sgm56*(p02(3,v-1)+p02(4,v-1))*lams2)/(h*c*Aco);
 W65=(Tp*sgm65p*(p02(5,v-1)+p02(6,v-1))*lamp+Ts2*sgm65*(p02(3,v-1)+p02(4,v-1))*lams2)/(h*c*Aco);
 
syms n2 n3 n6
eq1=-n2/t21+n3/t32+W12*(Ner-n2-n3)-W21*n2-2*C2*n2^2;%dN2/dt=0
eq2=-n3/t32+W13*(Ner-n2-n3)+R61*n6*(Ner-n2-n3)-R35*n3*(Nybp-n6)+C2*n2^2;%dN3/dt=0
eq3=-n6/t65+W56*(Nybp-n6)-W65*n6-R61*n6*(Ner-n2-n3)+R35*n3*(Nybp-n6);%dN6/dt=0

[N2 N3 N6]=solve(eq1,eq2,eq3,'n2','n3','n6');
N2=double(N2);
N3=double(N3);
N6=double(N6);
for i=1:3
    if N2(i)<Ner&&N2(i)>0&&N6(i)>0 %�������⣬ֻ��һ������0��������
        N2=N2(i);
        N3=N3(i);
        N6=N6(i);
        break
    end
end
%---nonparticipatory Yb3+ -ions--------------------------------------------
syms n6np
eq4=-n6np/t65+W56*(Nybnp-n6np)-W65*n6np;%dN6p/dt=0

[N6np]=solve(eq4,'n6np');
N6np=double(N6np);
%----------------ode45 ��� ���䷽��----------------------------------------
z02v=[L/m*(v-2) L/m*(v-1)];
p02v=p02(:,v-1);
[z2,y2]=ode45(@divide_double_pumped_fnp,z02v,p02v);

x2=length(z2(:,1));
Z2(v)=z2(x2); 
p2(:,v)=y2(x2,:);

p02(:,v)=p2(:,v);
end
%**************************************************************************
%------------------------p0(2,1)+eps2-----end------------------------------

%------------------------p0(6,1)+eps3-----begin----------------------------
%**************************************************************************
%----------------guess the initial values----------------------------------
% v=1; %z=0�����������
p03(5,1)=p0(5,1);           %pump+0
p03(6,1)=p0(6,1)+eps3;           %pump-0
p03(2,1)=p0(2,1);        %ps1-0
p03(4,1)=p0(4,1);             %ps2-0
p03(1,1)=R1er*p03(2,1);      %ps1+0
p03(3,1)=R1yb*p03(4,1);      %ps2+0
%-----------------p0(v)��L�ĵ�v�εĳ�ֵ,v=1��Z=0~1��һ������˵ĳ�ֵ��p0(2)�ǵڶ��εĳ�ֵ-----------------
%-----------------p(1)��v=1��Z=0~1��һ������˵ĳ�ֵ��p(2)����L��Z=0~1��һ���Ҷ˵ļ���ֵ��Ҳ�ǵڶ��εĳ�ֵ-----p(v+1)�ǵ�v�ε��Ҷ�ֵ----------------
Z3(1)=0; 
p3(:,1)=p03(:,1);
%----------------ode45��M������L���ȵĹ��˷ָ���M����һ��һ������һ��Z=Zv+1���ĳ�ֵ--------------

for v=2:m+1
%---------------��Z=z0v���� ���ܼ��������ֲ� n1 n2 n3 n5 n6------------------    
 W12=Ts1*sgm12*(p03(1,v-1)+p03(2,v-1))*lams1/(h*c*Aco);
 W21=Ts1*sgm21*(p03(1,v-1)+p03(2,v-1))*lams1/(h*c*Aco);
 W13=(Tp*sgm13p*(p03(5,v-1)+p03(6,v-1))*lamp+Ts2*sgm13*(p03(3,v-1)+p03(4,v-1))*lams2)/(h*c*Aco);
 W56=(Tp*sgm56p*(p03(5,v-1)+p03(6,v-1))*lamp+Ts2*sgm56*(p03(3,v-1)+p03(4,v-1))*lams2)/(h*c*Aco);
 W65=(Tp*sgm65p*(p03(5,v-1)+p03(6,v-1))*lamp+Ts2*sgm65*(p03(3,v-1)+p03(4,v-1))*lams2)/(h*c*Aco);
 
syms n2 n3 n6
eq1=-n2/t21+n3/t32+W12*(Ner-n2-n3)-W21*n2-2*C2*n2^2;%dN2/dt=0
eq2=-n3/t32+W13*(Ner-n2-n3)+R61*n6*(Ner-n2-n3)-R35*n3*(Nybp-n6)+C2*n2^2;%dN3/dt=0
eq3=-n6/t65+W56*(Nybp-n6)-W65*n6-R61*n6*(Ner-n2-n3)+R35*n3*(Nybp-n6);%dN6/dt=0

[N2 N3 N6]=solve(eq1,eq2,eq3,'n2','n3','n6');
N2=double(N2);
N3=double(N3);
N6=double(N6);
for i=1:3
    if N2(i)<Ner&&N2(i)>0&&N6(i)>0 %�������⣬ֻ��һ������0��������
        N2=N2(i);
        N3=N3(i);
        N6=N6(i);
        break
    end
end
%---nonparticipatory Yb3+ -ions--------------------------------------------
syms n6np
eq4=-n6np/t65+W56*(Nybnp-n6np)-W65*n6np;%dN6p/dt=0

[N6np]=solve(eq4,'n6np');
N6np=double(N6np);
%----------------ode45 ��� ���䷽��----------------------------------------
z03v=[L/m*(v-2) L/m*(v-1)];
p03v=p03(:,v-1);
[z3,y3]=ode45(@divide_double_pumped_fnp,z03v,p03v);

x3=length(z3(:,1));
Z3(v)=z3(x3); 
p3(:,v)=y3(x3,:);

p03(:,v)=p3(:,v);

end
%**************************************************************************
%------------------------p0(6,1)+eps3-----end------------------------------

%==========================================================================

%---D ΪF1��F2 F3  ��p0(2,1) p0(4,1), p0(6,1)�ĵ���-------------------------
D11=(p1(2,m+1)-p(2,m+1)-R2er*(p1(1,m+1)-p(1,m+1)))/eps1;
D12=(p2(2,m+1)-p(2,m+1)-R2er*(p2(1,m+1)-p(1,m+1)))/eps2;
D13=(p3(2,m+1)-p(2,m+1)-R2er*(p3(1,m+1)-p(1,m+1)))/eps3;
D21=(p1(4,m+1)-p(4,m+1)-R2yb*(p1(3,m+1)-p(3,m+1)))/eps1;
D22=(p2(4,m+1)-p(4,m+1)-R2yb*(p2(3,m+1)-p(3,m+1)))/eps2;
D23=(p3(4,m+1)-p(4,m+1)-R2yb*(p3(3,m+1)-p(3,m+1)))/eps3;
D31=(p1(6,m+1)-p(6,m+1))/eps1;
D32=(p2(6,m+1)-p(6,m+1))/eps2;
D33=(p3(6,m+1)-p(6,m+1))/eps3;

delta1=-(D12*D23*F3-D12*D33*F2-D13*D22*F3+D13*D32*F2+D22*D33*F1-D23*D32*F1)/(D11*D22*D33-D11*D23*D32-D12*D21*D33+D12*D23*D31+D13*D21*D32-D13*D22*D31);
delta2=(D11*D23*F3-D11*D33*F2-D13*D21*F3+D13*D31*F2+D21*D33*F1-D23*D31*F1)/(D11*D22*D33-D11*D23*D32-D12*D21*D33+D12*D23*D31+D13*D21*D32-D13*D22*D31);
delta3=-(D11*D22*F3-D11*D32*F2-D12*D21*F3+D12*D31*F2+D21*D32*F1-D22*D31*F1)/(D11*D22*D33-D11*D23*D32-D12*D21*D33+D12*D23*D31+D13*D21*D32-D13*D22*D31);

p0(2,1)=abs(p0(2,1)+delta1);             %ps1-0
p0(4,1)=abs(p0(4,1)+delta2);             %ps2-0
p0(6,1)=abs(p0(6,1)+delta3);             %pump-0
p0(1,1)=R1er*p0(2,1);      %ps1+0
p0(3,1)=R1yb*p0(4,1);      %ps2+0

end  %end for while
%==========================================================================

    ZZ_g(1,:,u,g)=Z;                   %z�ֳɵ�С�δ���һ��4ά����
    pp_g(:,:,u,g)=p;                   %p��z�仯������һ��4ά����

    Ti(g,u)=T;                         %core�¶ȴ���һ��������
    Pump(g,u)=p(5,1)+p(6,m+1);            %launched pump power����һ��2ά������
    Per(g,u)=p(1,m+1)+p(2,1);          %Er����������ȶ�ֵ����һ��2ά������
    Pyb(g,u)=p(3,m+1)+p(4,1);          %Yb����������ȶ�ֵ����һ��2ά������
    %---��ת�ó��У�������origin��ϣ���sigmaplot��ͼ---------------------------


%--------------------------------------------------------------------------

end  %end for p56

subplot(2,4,8);
 switch g
    case 1
        
        plot(Pump(g,:),Per(g,:),'r-s');hold on;%��Er���źŹ��������ͼ
        plot(Pump(g,:),Pyb(g,:),'b-s');hold on;%��Yb���źŹ��������ͼ
     
     case 2
      
        plot(Pump(g,:),Per(g,:),'r--o');hold on;%��Er���źŹ��������ͼ
        plot(Pump(g,:),Pyb(g,:),'b--o');hold on;%��Yb���źŹ��������ͼ

    case 3
        
        plot(Pump(g,:),Per(g,:),'r:^');hold on;%��Er���źŹ��������ͼ
        plot(Pump(g,:),Pyb(g,:),'b:^');hold on;%��Yb���źŹ��������ͼ
    
   
    otherwise
    break;
 end

end %end for T
    title('Thermal Effect to Er,Yb co-doped Fiber Laser Source');
    xlabel('Launched pump power/W');
    ylabel('Output power/W');
    legend('1550nm','1064nm 273K','1550nm','1064nm 323K','1550nm','1064nm 373K','Location','NorthWest');
%------------------------------------------------------------------------------------------------------------------   
    Pump_rot=rot90(Pump);
        Per_rot=rot90(Per);
            Pyb_rot=rot90(Pyb);
%------------------------------------------------------------------------------------------------------------------   
            
subplot(2,4,7);
plot(Z,p(1,:),'rs',Z,p(2,:),'r*',Z,p(3,:),'b^',Z,p(4,:),'b+',Z,p(5,:),'g.-',Z,p(6,:),'g.:');hold on;
title('Er,Yb co-doped fiber Laser source');
xlabel('Position of the fiber/m');
ylabel('Power/W');
legend('Er 1550nm +','Er 1550nm -','Yb 1064nm +','Yb 1064nm- ','Pump 975nm +','Pump 975nm -','Location','North');

toc
%===================================================================================================================

switch NNNN
    case 1
        save 2013R001
    case 2
        save 2013R0005
    %case 3
    %    save 2013Ner7.6
    %case 4
     %   save 2013Ner6.6
        otherwise
    break;
end 
end


