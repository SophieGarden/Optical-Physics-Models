
clear all;
clc;
close all;

addpath('C:/Users/Wendy/documents/Research/Hg/HgData/HgCellwithVaccum/248nmPumped/11-04-2015/237C')

Co = zeros(3,5);
for j=1:2;%j=1 for 335; 2 for 485
    figure
    for i= 1:5;
        if j==1
            if i==1
                M = csvread('335nm-1.csv',23,0);
            elseif i==2
                M = csvread('335nm-2.csv',23,0);
            elseif i==3
                M = csvread('335nm-3.csv',23,0);
            elseif i==4
                M = csvread('335nm-4.csv',23,0);
            elseif i==5
                M = csvread('335nm-5.csv',23,0);
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif j==2
            if i==1
                M = csvread('485nm-1.csv',23,0);
            elseif i==2
                M = csvread('485nm-2.csv',23,0);
            elseif i==3
                M = csvread('485nm-3.csv',23,0);
            elseif i==4
                M = csvread('485nm-4.csv',23,0);
            elseif i==5
                M = csvread('485nm-5.csv',23,0);
            end
            %
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
        end
        x = M(:,1)*1e6; %unit: micro second
        y = - M(:,3);
        
        plot (x,y,'g'); hold on;
        if j==1
            x1 = x((x>=2) );  %x1 for fit  %x>=5 for 335nm x>=10 for 485nm
            y1 = y((x>=2) );  %y1 for fit
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0 0.001 0],...
                'Upper',[1 0.1 0.1],...
                'StartPoint',[0.0046 0.0171 0.0041]);
            
        elseif j==2
            x1 = x((x>=55) );  %x1 for fit  When the decay begins: x>=15 for 485nm
            y1 = y((x>=55) );  %y1 for fit
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0 0.001 0],...
                'Upper',[1 0.1 0.1],...
                'StartPoint',[0.0147 0.0110 0.0044]);
            
        end
        plot (x1,y1,'b');hold on;
        
        %     fo = fitoptions('Method','NonlinearLeastSquares',...
        %         'Lower',[0 0.001 0],...
        %         'Upper',[1 100 0.1],...
        %         'StartPoint',[2.37e-5 30 0.0064]);
        
        myfittype = fittype('a*exp(-x*k)+c',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a','k','c'},'options',fo);
        
        myfit = fit(x1,y1,myfittype);
        Co(:,i) = coeffvalues(myfit);
        plot(myfit,x1,y1);hold on;
    end
    if j==1
        Co2(1:3,:) = Co(:,:); %decay rate unit: 1e6 sec-1
        Co2(4,:) = 1./Co(2,:); %lifetime unit:us
        Co3(:,:) = Co2(:,:);
        Co3(:,6) = mean(Co2(:,:),2);
        Co3(:,7) = std(Co2,1,2);
        Co4(1:4,:) = Co3;
    elseif j==2
        Co2(1:3,:) = Co(:,:); %decay rate unit: 1e6 sec-1
        Co2(4,:) = 1./Co(2,:); %lifetime unit:us
        Co3(:,:) = Co2(:,:);
        Co3(:,6) = mean(Co2(:,:),2);
        Co3(:,7) = std(Co2,1,2);
        Co4(5:7,:) = Co3;
    end
end


