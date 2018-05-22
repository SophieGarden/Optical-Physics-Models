function dn=thulium330(t,n)
%或者全体都不代表光子数密度
% n1代表基态粒子数密度；n2代表激光高能态粒子数密度；n代表反转粒子数密度；fai代表谐振腔内光子数密度
%以下所有长度单位为米

sigmaas=0.56e-25;
sigmaes=5.96e-25;%受激发射截面
sigmaap=1.52e-25;%受激吸收截面
sigmaep=0.11e-25;
gamal=1+sigmaas/sigmaes;%gamal=1+fll/ful
gamap=1+sigmaep/sigmaap;%gamap=1+fup/flp
r=(13e-6)/2;%纤芯半径
A=pi*r^2;%纤芯面积
l1=0.3;%增益介质的长度
%deltas=2.3e-3;
%deltap=1.2e-2;
ntot=0.4e26; %掺杂浓度
c0=3e8;
nindex=1.46;%纤芯折射率
c=c0/nindex;
lanbdap=1.565e-6;%泵浦波长
lanbdas=1.86e-6;%激光波长
h=6.626e-34;
tao21=334e-6;%激光高能级寿命
l=0.3;%谐振腔的光学长度
R1=1;%左腔镜对激光的反射率
R2=0.04;%右腔镜对激光的反射率
alfa=2.3e-3;%腔损耗（衍射吸收等）
deltac=-log(R1*R2)+2*alfa*l;
taoc=2*l/(c*deltac);

deltaN=1/(taoc*c0*sigmaes);%阈值时的反转粒子数

deltat=1e-6;
P0=6;%泵浦光幅值
%梯形脉冲  
t1=0.10e-6;
t2=0.7e-6;
t3=0.8e-6;
t4=t1+deltat;
t5=t2+deltat;
t6=t3+deltat;
t7=t1+2*deltat;
t8=t2+2*deltat;
t9=t3+2*deltat;
t10=t1+3*deltat;
t11=t2+3*deltat;
t12=t3+3*deltat;
t13=t1+4*deltat;
t14=t2+4*deltat;
t15=t3+4*deltat;
t16=t1+5*deltat;
t17=t2+5*deltat;
t18=t3+5*deltat;
t19=t1+6*deltat;
t20=t2+6*deltat;
t21=t3+6*deltat;
t22=t1+7*deltat;
t23=t2+7*deltat;
t24=t3+7*deltat;
t25=t1+8*deltat;
t26=t2+8*deltat;
t27=t3+8*deltat;
t28=t1+9*deltat;
t29=t2+9*deltat;
t30=t3+9*deltat;
t31=t1+10*deltat;
t32=t2+10*deltat;
t33=t3+10*deltat;
t34=t1+11*deltat;
t35=t2+11*deltat;
t36=t3+11*deltat;
t37=t1+12*deltat;
t38=t2+12*deltat;
t39=t3+12*deltat;
t40=t1+13*deltat;
t41=t2+13*deltat;
t42=t3+13*deltat;
t43=t1+14*deltat;
t44=t2+14*deltat;
t45=t3+14*deltat;
t46=t1+15*deltat;
t47=t2+15*deltat;
t48=t3+15*deltat;
t49=t1+16*deltat;
t50=t2+16*deltat;
t51=t3+16*deltat;
t52=t1+17*deltat;
t53=t2+17*deltat;
t54=t3+17*deltat;
t55=t1+18*deltat;
t56=t2+18*deltat;
t57=t3+18*deltat;
t58=t1+19*deltat;
t59=t2+19*deltat;
t60=t3+19*deltat;

%J=6e-6;%泵浦脉冲能量，单位焦耳
%P0=J/((t2+t3)/2-t1/2);%P为泵浦功率
if t<t1
    P=P0*t/t1;
else
    if t>=t1&t<=t2
        P=P0;
    else
            if t>t2&t<t3
                P=(t-t3)*P0/(t2-t3);
            else 
                if t>=t3&t<=deltat
                    P=0;
            else
if t>deltat&t<t4
    P=P0*(t-deltat)/(t4-deltat);
else
    if t>=t4&t<=t5
        P=P0;
    else
            if t>t5&t<t6
                P=(t-t6)*P0/(t5-t6);
            else 
                if t>=t6&t<=2*deltat
                P=0;
                
                else
if t>2*deltat&t<t7
    P=P0*(t-2*deltat)/(t7-2*deltat);
else
    if t>=t7&t<=t8
        P=P0;
    else
            if t>t8&t<t9
                P=(t-t9)*P0/(t8-t9);
            else 
                if t>=t9&t<=3*deltat
                P=0;
                else
if t>3*deltat&t<t10
    P=P0*(t-3*deltat)/(t10-3*deltat);
else
    if t>=t10&t<=t11
        P=P0;
    else
            if t>t11&t<t12
                P=(t-t12)*P0/(t11-t12);
            else 
                if t>=t12&t<=4*deltat
                P=0;
                else
if t>4*deltat&t<t13
    P=P0*(t-4*deltat)/(t13-4*deltat);
else
    if t>=t13&t<=t14
        P=P0;
    else
            if t>t14&t<t15
                P=(t-t15)*P0/(t14-t15);
            else 
                if t>=t15&t<=5*deltat
                P=0;
                else
if t>5*deltat&t<t16
    P=P0*(t-5*deltat)/(t16-5*deltat);
else
    if t>=t16&t<=t17
        P=P0;
    else
            if t>t17&t<t18
                P=(t-t18)*P0/(t17-t18);
            else 
                if t>=t18&t<=6*deltat
                P=0;
                else
if t>6*deltat&t<t19
    P=P0*(t-6*deltat)/(t19-6*deltat);
else
    if t>=t19&t<=t20
        P=P0;
    else
            if t>t20&t<t21
                P=(t-t21)*P0/(t20-t21);
            else 
                if t>=t21&t<=7*deltat
                P=0;
                else
if t>7*deltat&t<t22
    P=P0*(t-7*deltat)/(t22-7*deltat);
else
    if t>=t22&t<=t23
        P=P0;
    else
            if t>t23&t<t24
                P=(t-t24)*P0/(t23-t24);
            else 
                if t>=t24&t<=8*deltat
                P=0;
                else
if t>8*deltat&t<t25
    P=P0*(t-8*deltat)/(t25-8*deltat);
else
    if t>=t25&t<=t26
        P=P0;
    else
            if t>t26&t<t27
                P=(t-t27)*P0/(t26-t27);
            else 
                if t>=t27&t<=9*deltat
                P=0;
                else
if t>9*deltat&t<t28
    P=P0*(t-9*deltat)/(t28-9*deltat);
else
    if t>=t28&t<=t29
        P=P0;
    else
            if t>t29&t<t30
                P=(t-t30)*P0/(t29-t30);
            else
                if t>=t30&t<=10*deltat
                P=0;
                else
                
if t>10*deltat&t<t31
    P=P0*(t-10*deltat)/(t31-10*deltat);
else
    if t>=t31&t<=t32
        P=P0;
    else
            if t>t32&t<t33
                P=(t-t33)*P0/(t32-t33);
            else 
                if t>=t33&t<=11*deltat
                    P=0;
            else
if t>11*deltat&t<t34
    P=P0*(t-11*deltat)/(t34-11*deltat);
else
    if t>=t34&t<=t35
        P=P0;
    else
            if t>t35&t<t36
                P=(t-t36)*P0/(t35-t36);
            else 
                if t>=t36&t<=12*deltat
                P=0;
                
                else
if t>12*deltat&t<t37
    P=P0*(t-12*deltat)/(t37-12*deltat);
else
    if t>=t37&t<=t38
        P=P0;
    else
            if t>t38&t<t39
                P=(t-t39)*P0/(t38-t39);
            else 
                if t>=t39&t<=13*deltat
                P=0;
                else
if t>13*deltat&t<t40
    P=P0*(t-13*deltat)/(t40-13*deltat);
else
    if t>=t40&t<=t41
        P=P0;
    else
            if t>t41&t<t42
                P=(t-t42)*P0/(t41-t42);
            else 
                if t>=t42&t<=14*deltat
                P=0;
                else
if t>14*deltat&t<t43
    P=P0*(t-14*deltat)/(t43-14*deltat);
else
    if t>=t43&t<=t44
        P=P0;
    else
            if t>t44&t<t45
                P=(t-t45)*P0/(t44-t45);
            else 
                if t>=t45&t<=15*deltat
                P=0;
                else
if t>15*deltat&t<t46
    P=P0*(t-15*deltat)/(t46-15*deltat);
else
    if t>=t46&t<=t47
        P=P0;
    else
            if t>t47&t<t48
                P=(t-t48)*P0/(t47-t48);
            else 
                if t>=t48&t<=16*deltat
                P=0;
                else
if t>16*deltat&t<t49
    P=P0*(t-16*deltat)/(t49-16*deltat);
else
    if t>=t49&t<=t50
        P=P0;
    else
            if t>t50&t<t51
                P=(t-t51)*P0/(t50-t51);
            else 
                if t>=t51&t<=17*deltat
                P=0;
                else
if t>17*deltat&t<t52
    P=P0*(t-17*deltat)/(t52-17*deltat);
else
    if t>=t52&t<=t53
        P=P0;
    else
            if t>t53&t<t54
                P=(t-t54)*P0/(t53-t54);
            else 
                if t>=t54&t<=18*deltat
                P=0;
                else
if t>18*deltat&t<t55
    P=P0*(t-18*deltat)/(t55-18*deltat);
else
    if t>=t55&t<=t56
        P=P0;
    else
            if t>t56&t<t57
                P=(t-t57)*P0/(t56-t57);
            else 
                if t>=t57&t<=19*deltat
                P=0;
                else
if t>19*deltat&t<t58
    P=P0*(t-19*deltat)/(t58-19*deltat);
else
    if t>=t58&t<=t59
        P=P0;
    else
            if t>t59&t<t60
                P=(t-t60)*P0/(t59-t60);
            else 
                P=0;
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
               end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end
                end
            end
    end
end


flp=1;
fup=sigmaep/sigmaap;
%Rpp=0.04;

%Rp0=P/(A*h*c0/lanbdap);
%Rp=Rp0*(1-exp(-((flp*(ntot-n(1))-fup*n(1))*sigmaap*l1)))* ...
%    (1+Rpp*exp(-((flp*(ntot-n(1))-fup*n(1))*sigmaap*l1)));
%dn(1)=Rp/l1-c*sigmaes*(gamal*n(1)-(gamal-1)*ntot)*n(2)-n(1)/tao21; %N2能级粒子数密度变化率
%dn(2)=c*sigmaes*(l1/l)*(gamal*n(1)-(gamal-1)*ntot)*n(2)-n(2)/taoc+0.3e-3*n(1)/tao21;%激光谐振腔内光子数密度变化率fai

dn(1)=sigmaap*P*(ntot-gamap*n(1))/(A*h*c0/lanbdap)-c*sigmaes*(gamal*n(1)-(gamal-1)*ntot)*n(2)-n(1)/tao21; %N2能级粒子数密度变化率
dn(2)=c*sigmaes*(l1/l)*(gamal*n(1)-(gamal-1)*ntot)*n(2)-n(2)/taoc+0.3e-3*n(1)/tao21;%激光谐振腔内光子数密度变化率fai


dn=[dn(1);dn(2)];
