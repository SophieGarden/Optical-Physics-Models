n1=1.550   %��о������========%����������
n2=1.536   %�ڰ���������
a=25e-6 %��о�뾶
lamda1=1e-6  % 1��m����
lamda2=1.55e-6 % 1.55��m����

k1=2*pi/lamda1
k2=2*pi/lamda2

delta=(n1^2-n2^2)/2./n1^2 %��һ��Ƶ��

V1=k1*a*(n1^2-n2^2)^0.5 
V2=k2*a*(n1^2-n2^2)^0.5 
%%%%%%%%%%
v=0
ev=2 %v=0%ev=2%v<>0%ev=1
m=v
n=1   
P=m+n+1 % LPmn

bata1=n1*k1*(1-2*(2*delta)^0.5*P/a/n1/k1)^0.5
bata2=n1*k2*(1-2*(2*delta)^0.5*P/a/n1/k2)^0.5

kapa1=((n1*k1)^2-bata1^2)^0.5
kapa2=((n1*k2)^2-bata2^2)^0.5

gama1=(bata1^2-(n2*k1)^2)^0.5
gama2=(bata2^2-(n2*k2)^2)^0.5
%%%%��ֵ����%%%

%K(v-1)=BESSELK(gama*a,v-1)
%K(v+1)=BESSELK(gama*a,v+1)



R=0.004:0.00001:0.005

Alfa1=2*a*kapa1^2*exp(2*a*gama1)*exp(-2/3*gama1^3/bata1^2.*R)/ev/(pi*gama1)^0.5/V1^2./R.^0.5
Alfa2=2*a*kapa2^2*exp(2*a*gama2)*exp(-2/3*gama2^3/bata2^2.*R)/ev/(pi*gama2)^0.5/V2^2./R.^0.5

plot(R,Alfa1,'r',R,Alfa2,'b')
xlabel('R');ylabel('2��')




