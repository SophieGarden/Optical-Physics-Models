%%fit for temperature dependence of coefficients
close all;
clear;
clc;
x0=[107,138,152,162,173,182,195,204,218,229,237,266,284,298,309,319,327,335,341,348] + 273.15;
y0=[1000,80,30,25,15,3,1.20000000000000,1,0.500000000000000,0.300000000000000,0.200000000000000,...
    0.0500000000000000,0.0250000000000000,0.0150000000000000,0.00900000000000000,...
    0.00600000000000000,0.00600000000000000,0.00400000000000000,0.00350000000000000,...
    0.00300000000000000];
x1=transpose(x0);
y1=transpose(y0);
% y1 = log(y1);
w=1./y1.^2;

    fo = fitoptions('Method','NonlinearLeastSquares',... 
        ...%         'Lower',[0 0.001 0],...
...%         'Upper',[10 0.2 0.1],...
        'StartPoint',[1.00000000000000e+16,0.0787457109666162,1.06590151967650]);
fo.weight= w;
 
plot (x1,y1,'b');hold on;

myfittype = fittype('a*exp(-x*k)+c',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','k','c'},'options',fo);

% myfittype = fittype('poly2');
myfit = fit(x1,y1,myfittype);
Co = coeffvalues(myfit);
plot(myfit,x1,y1);hold on;
plot(x1,myfit(x1),'ro');
% exp decay
 y2=myfit(x1)-y1;
%poly2
%y2=exp(myfit(x1))-y1;

y3=y2./y1;

%% fIt for high T
%fit for temperature dependence of coefficients
close all;
clear;
clc;
x0=[107,138,152,162,173,182,195,204,218,229,237,266,284,298,309,319,327,335,341,348] + 273.15;
y0=[1000,80,30,25,15,3,1.20000000000000,1,0.500000000000000,0.300000000000000,0.200000000000000,...
    0.0500000000000000,0.0250000000000000,0.0150000000000000,0.00900000000000000,...
    0.00600000000000000,0.00600000000000000,0.00400000000000000,0.00350000000000000,...
    0.00300000000000000];
x1=transpose(x0);
y1=transpose(y0);
x1 = x1(1:20);
y1=y1(1:20);
w=x1;

fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1.00000000000000e+16,0.0787451125340310,0.820225530518243]);
fo.weight= w;
 
plot (x1,y1,'b');hold on;
myfittype = fittype('a*exp(-x*k)+c',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','k','c'},'options',fo);

% myfittype = fittype('poly2');
myfit = fit(x1,y1,myfittype);
Co = coeffvalues(myfit);
plot(myfit,x1,y1);hold on;
plot(x1,myfit(x1),'ro');

y2=myfit(x1)-y1;

y3=y2./y1;
