
clear all;
clc;
close all;
%addpath('HgCellwithVaccum/266nmPumped/06-25-2015/173C')
% addpath('HgCellwithVaccum/266nmPumped/07-22-2015/218C')
order = 16;
    if order ==1
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/107C')
    elseif order ==2
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/138C')
    elseif order ==3
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/152C')
    elseif order ==4
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/162C')
    elseif order ==5
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/173C')
    elseif order ==6
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/182C')
    elseif order ==7
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/195C')
    elseif order ==8
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/204C')
    elseif order ==9
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/218C')
    elseif order ==10
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/229C')
    elseif order ==11
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/237C')
    elseif order ==12
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/266C')
    elseif order ==13
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/284C')
    elseif order ==14
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/298C')
    elseif order ==15
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/309C')
    elseif order ==16
        addpath('C:\Users\TotoroMuscles\Documents\Research\Hg\HgData_Processing\HgCellwithVaccum\266nmPumped\finaldata\319C')
    elseif order ==17
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/327C')
    elseif order ==18
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/335C')
    elseif order ==19
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/341C')
    elseif order ==20
        addpath('C:/Users/Wendy/documents/Research/Hg/HgData_Processing/HgCellwithVaccum/266nmPumped/finaldata/348C')
    end
 Co = zeros(3,5);
    j=2;%j=1 for 335; 2 for 485
    for i= 1:5;
        if j==1
            if i==1
                M = csvread('335nm-200uj-1.csv',23,0);
            elseif i==2
                M = csvread('335nm-200uj-2.csv',23,0);
            elseif i==3
                M = csvread('335nm-200uj-3.csv',23,0);
            elseif i==4
                M = csvread('335nm-200uj-4.csv',23,0);
            elseif i==5
                M = csvread('335nm-200uj-5.csv',23,0);
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif j==2
            if i==1
                M = csvread('485nm-20uj-1.csv',23,0);
            elseif i==2
                M = csvread('485nm-20uj-2.csv',23,0);
            elseif i==3
                M = csvread('485nm-20uj-3.csv',23,0);
            elseif i==4
                M = csvread('485nm-20uj-4.csv',23,0);
            elseif i==5
                M = csvread('485nm-20uj-5.csv',23,0);
            end
            %
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            %         if i==1
            %         M = csvread('20uj-485-1.csv',23,0);
            %     elseif i==2
            %         M = csvread('20uj-485-2.csv',23,0);
            %     elseif i==3
            %         M = csvread('20uj-485-3.csv',23,0);
            %     elseif i==4
            %         M = csvread('20uj-485-4.csv',23,0);
            %     elseif i==5
            %         M = csvread('20uj-485-5.csv',23,0);
            %         end
        end
        x = M(:,1)*1e6; %unit: micro second
        %     y = - M(:,2);
        span=101;
        y = smooth(- M(:,2),span) ;
        plot (x,y,'g'); hold on;
        if j==1
            x1 = x((x>=2) );  %x1 for fit  %x>=5 for 335nm x>=10 for 485nm
            y1 = y((x>=2) );  %y1 for fit
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0 0.001 0],...
                'Upper',[1 0.1 0.1],...
                'StartPoint',[0.03 0.055 0.0101]);
            
        elseif j==2
            b=15;
            x1 = x((x>=b) );  %x1 for fit  When the decay begins: x>=15 for 485nm
            y1 = y((x>=b) );  %y1 for fit
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0 0.001 0],...
                'Upper',[1 0.1 0.1],...
                'StartPoint',[0.003 0.055 0.01]);
            
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
    
    Co2(1:3,:) = Co(:,:); %decay rate unit: 1e6 sec-1
    Co2(4,:) = 1./Co(2,:); %lifetime unit:us
    Co3(:,:) = Co2(:,:);
    Co3(:,6) = mean(Co2(:,:),2);
    
    Co3(:,7) = std(Co2,1,2);
       Co3(:,8) = b;
 
