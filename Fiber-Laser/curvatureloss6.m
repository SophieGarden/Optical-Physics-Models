%bend loss varies with NA
function y=curvatureloss6
n1=1.55;
n2=1.536;% 没有考虑到n与波长lamd有关
a=25e-6;%米
m=0;n=1;
%R=8e-3;
ev=2;
for c=1:4
    na=0.2+(c-1)*0.01;
    n2=sqrt(n1^2-na^2);
%NA(c)=na;
lamd=[1.064e-6,1.55e-6];

for j=1:2
    i=1;
for r=(4.65:0.005:5.2)*1e-3-(c-1)*0.4e-3;
    
k=2*3.1415926/lamd(j);
V=k*a*sqrt(n1^2-n2^2);
delta=(n1^2-n2^2)/(2*n1^2);
%bat=k*n1*sqrt(1-2*sqrt(2*delta)*(m+n+1)/(k*n1*a));
bat=n2*k*(0.95*delta+1);   
gama=sqrt((bat^2-n2^2*k^2));
ki=sqrt(n1^2*k^2-bat^2);
loss=2*a*ki^2*exp(2*gama*a-2/3*gama^3/bat^2*r)/(ev*sqrt(3.1415926*gama*r)*V^2);
LOSS(j,i)=loss;
R(i)=r;
i=i+1;
end  
end
switch c
    case 1
        subplot(2,2,1);
        plot(R,LOSS(1,:),'-g');hold on;
        plot(R,LOSS(2,:),'-r');hold on;
        xlabel('R');ylabel('2α');
        legend('1.064μm','1.55μm');
        title('NA=0.2');
    case 2
        subplot(2,2,2);
        plot(R,LOSS(1,:),'-g');hold on;
        plot(R,LOSS(2,:),'-r');hold on;
        xlabel('R');ylabel('2α');
        legend('1.064μm','1.55μm');
        title('NA=0.21');
    case 3
        subplot(2,2,3);
        plot(R,LOSS(1,:),'-g');hold on;
        plot(R,LOSS(2,:),'-r');hold on;
        xlabel('R');ylabel('2α');
        legend('1.064μm','1.55μm');
        title('NA=0.22');
    case 4
        subplot(2,2,4);
        plot(R,LOSS(1,:),'-g');hold on;
        plot(R,LOSS(2,:),'-r');hold on;
        xlabel('R');ylabel('2α');
        legend('1.064μm','1.55μm');
        title('NA=0.23');
    otherwise
        break;
end

end

end