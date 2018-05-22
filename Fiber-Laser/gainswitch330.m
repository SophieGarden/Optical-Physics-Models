function gainswitch330   
%根据沈老师的实验参数进行模拟

sigmaas=0.56e-25;
sigmaes=5.96e-25;%受激发射截面
sigmaap=1.52e-25;%受激吸收截面
sigmaep=0.11e-25;
gamal=1+sigmaas/sigmaes;%gamal=1+fll/ful
gamap=1+sigmaep/sigmaap;%gamap=1+fup/flp
c0=3e8;
nindex=1.46;%纤芯折射率
c=c0/nindex;
lanbdap=1.565e-6;%泵浦波长
lanbdas=1.86e-6;%激光波长
h=6.626068e-34;
r=(13e-6)/2;%纤芯半径
A=pi*r^2;%纤芯面积
l1=0.3;%增益介质的长度
l=0.3;%谐振腔的光学长度

R1=1;%左腔镜对激光的反射率
R2=0.04;%右腔镜对激光的反射率
alfa=2.3e-3;%腔损耗（衍射吸收等）
ntot=0.4e26; %掺杂浓度/总的粒子数密度
H=[0,20e-6];n0=[ntot*1e-11,1e3];[t,n]=ode15s('thulium330',H,n0);

deltac=-log(R1*R2)+2*alfa*l;
taoc=2*l/(c*deltac);
%Ps(:,2)=n(:,2)*A*l1*(h*c0/lanbdas)/taoc;
%Ps(:,2)=A*c*n(:,2)*(h*c/lanbdas);

Pout(:,2)=-A*c*n(:,2)*(h*c0/lanbdas)*log(R2)/2;
%Pout(:,2)=(1/2)*c*n(:,2)*h*(c0/lanbdas)*A*(l1/l)*log(1/R2);
Pmax=max(Pout(:,2))
fito=trapz(t, Pout(:,2))%total photon numbers inside laser cavity-one pulse
FWHM=(fito/Pmax)*1e9


plot(t,Pout(:,2),'g-')
hold on

%泵浦梯形脉冲
deltat=1e-6;%一个周期的时间
P0=6;
for i=0:1:19  %脉冲数：按20个脉冲计算
t1=0.10e-6+i*deltat;
t2=0.7e-6+i*deltat;
t3=0.8e-6+i*deltat;
%J=6e-6;%泵浦脉冲能量，单位焦耳
%P0=J/((t2+t3)/2-t1/2);%P为泵浦功率

t=i*deltat:1e-9:t1;
P=P0*(t-i*deltat)/(t1-i*deltat); %pump
plot(t,P,'y-')
hold on
t=t1:1e-8:t2;
P=P0;
plot(t,P,'y-')
hold on
t=t2:1e-9:t3;
P=(t-t3)*P0/(t2-t3);
plot(t,P,'y-')
hold on
t=t3:1e-8:(i+1)*deltat;
P=0;  
plot(t,P,'y-')
hold on;
ylabel('Pump (W)');
end






