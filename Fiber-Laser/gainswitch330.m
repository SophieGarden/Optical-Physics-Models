function gainswitch330   
%��������ʦ��ʵ���������ģ��

sigmaas=0.56e-25;
sigmaes=5.96e-25;%�ܼ��������
sigmaap=1.52e-25;%�ܼ����ս���
sigmaep=0.11e-25;
gamal=1+sigmaas/sigmaes;%gamal=1+fll/ful
gamap=1+sigmaep/sigmaap;%gamap=1+fup/flp
c0=3e8;
nindex=1.46;%��о������
c=c0/nindex;
lanbdap=1.565e-6;%���ֲ���
lanbdas=1.86e-6;%���Ⲩ��
h=6.626068e-34;
r=(13e-6)/2;%��о�뾶
A=pi*r^2;%��о���
l1=0.3;%������ʵĳ���
l=0.3;%г��ǻ�Ĺ�ѧ����

R1=1;%��ǻ���Լ���ķ�����
R2=0.04;%��ǻ���Լ���ķ�����
alfa=2.3e-3;%ǻ��ģ��������յȣ�
ntot=0.4e26; %����Ũ��/�ܵ��������ܶ�
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

%������������
deltat=1e-6;%һ�����ڵ�ʱ��
P0=6;
for i=0:1:19  %����������20���������
t1=0.10e-6+i*deltat;
t2=0.7e-6+i*deltat;
t3=0.8e-6+i*deltat;
%J=6e-6;%����������������λ����
%P0=J/((t2+t3)/2-t1/2);%PΪ���ֹ���

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






