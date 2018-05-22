% determine the vibrational constants We, WeXe; subscript u for upper
% electronic state, l for lower electronic state


% ATTENTION!
%!!!!!!!!!!!!!!!!! m n is the order of upper state!!!!!!!!!!!!!!!!!!!!!!!!!!!!


clear all;
clc;
% tic

Vl = 0:1:20; % V1 = number of peaks -1
peakno = length(Vl);
x1 = zeros(1,peakno);
x2 = zeros(1,peakno);
Y1 = zeros(1,peakno);
Y2 = zeros(1,peakno);

kk = 3; % larger than 1
mm = 3;

a = 8888;

coefficient1 = a(ones(mm+1,3));

diff1 =a(ones(mm+1));

kkkk=kk+1-max(-kk,-mm);
A = a(ones(kkkk,1));

Weuf = A;
WeXeuf = A;
Vef = A;

Wel = a(ones(kkkk,1));
WeXel = a(ones(kkkk,1));
Vel = a(ones(kkkk,1));

diff22 = A;
diff33 = A;
RMSD = A;

K = A;
El = a(ones(kkkk,15));
Eu =  a(ones(kkkk,15));

DE1 = zeros(1,peakno);
DE2 = zeros(1,peakno);

B = a(ones(kkkk,peakno));
DDE1 =  a(ones(kkkk,peakno));
DDE2 =  a(ones(kkkk,peakno));
DDE10 =  a(ones(kkkk,peakno));

DE1s = zeros(kkkk,peakno+5);
DE2s = zeros(kkkk,peakno+5);
DE3s = zeros(kkkk,peakno+5);
DE4s = zeros(kkkk,peakno+5);

Total=a(ones((mm+1)*kkkk,9));

syms Weu1 WeXeu1 Weu2 WeXeu2 Ve1 Ve2

load('y14200cmGB1V4.mat')
load('y14200cmGB2V4.mat')
load('y21550cmGB1V4.mat')
load('y21550cmGB2V4.mat')

load('y16720.mat'); % y1 fit peaks :data difference very big...not as good as find peaks values
load('y17702.mat'); % y2
load('y21126.mat');
load('y23483.mat');





for j = 1:peakno
       

% ---------14200 cm-1--------

%     Y1(j) = y14200cmGB1V4(j); %33
%     Y1(j) = y14200cmGB2V4(j); %17

%-----21550cm-1----------
% Y1(j) = y21550cmGB1V4(j); %10
% Y1(j) = y21550cmGB2V4(j); %13

%-----16720cm-1----------
% Y1(j) = y16720(j); %21

%-----17702cm-1----------
%Y1(j) = y17702(j); %21

%-----21126 cm-1---------
  Y1(j) = y21126(j); %25

end

ii=1;
for i = 1:5
    %===============================================================
    if i == 1
        %X1 sigma+, ground state, singlet
        Vel(i) = 0;
        Wel(i) = 215.15; %ground state constants
        WeXel(i)= 0.411;
    elseif i == 2
        %D1 Pia, singlet
        Vel(i) = 20689.6;
        Wel(i) = 199.2;
        WeXel(i)= 0.44;
    elseif i == 3
        %C1 Sigma,singlet
        Vel (i) = 12401.5062;
        Wel(i) = 187.4554;
        WeXel(i) = 0.4353; 
    elseif i == 4
        %  B1 Pia, singlet
        Vel(i) = 9917.3;
        Wel (i)= 192.21;
        WeXel (i) = 0.463;
    elseif i == 5
        %  some Omega, singlet
        Vel(i) = 24107.15;
        Wel(i) = 196.9;
        WeXel(i) = 0.4;
    
        
    end
    %===============================================================

    for m = 0:1:mm % m is the order of transition of the UPPER electronic state of the First group
        for j = 1:peakno
            x1(j) = Vl(j)+m; % x1 is v' of the First group
        end
        p1 = polyfit(x1,Y1,2);
                coefficient1(m+1,:) = p1;
        
        for k = -kk:1:min(kk,m)   % k is difference of the order of transition vibrational states between upper and lower electronic transitions
            % for m=1,deltaV1 min=-1, n=0, deltaV2 min=0
            
            deltaV1 = k;     % Vu = Vl+deltaV; deltaV is constant for the three groups
            %E(vl) = Ve+(Weu*(Vu+0.5)-WeXeu*(Vu+0.5)^2)-(Wel*(Vl+0.5)-WeXel*(Vl+0.5)^2); %transition emitted energy
            deltaV2 = deltaV1 + 1;
            deltaV3 = deltaV1 + 2;
            deltaV4 = deltaV1 - 1;
            
            
            % group 1
            A1 = [0 -1 0; 1 -1 0; 0.5 -0.25 1];
            b1 = [p1(1)-WeXel(i); p1(2)+Wel(i)-WeXel(i)*(1-2*deltaV1); p1(3)+Wel(i)*(0.5-deltaV1)-WeXel(i)*(0.5-deltaV1)^2];
            X1 = A1\b1;
            
            % final results of these vibrational constants corresponding to each m, n and k
%             Weuf(ii) = X1(1); % final results of these vibrational constants
%             WeXeuf(ii) = X1(2);
%             Vef(ii) = X1(3);
            
             % final results of these vibrational constants corresponding to each m, n and k
            Weuf(ii) = round(X1(1)*1000)/1000; % final results of these vibrational constants
            WeXeuf(ii) = round(X1(2)*1000)/1000;
            Vef(ii) = round(X1(3)*1000)/1000;
%         
        K(ii) = k; %record value of k
            %===============================================================================================================
            for s=0:1:peakno+m+k+25
                v = s;
                % ground state calculated energy value % El(1,2,3,...)
                % corresponding to V=0,1,2,...
                El(s+1,ii) = Vel(i)+Wel(i)*(v+0.5)-WeXel(i)*(v+0.5)^2;
                Eu(s+1,ii) = Vel(i)+Vef(ii)+Weuf(ii)*(v+0.5)-WeXeuf(ii)*(v+0.5)^2;
            end
            
            for ss = 0:1:peakno-1
                % for first group position of line in cm-1
                DE1(ii,ss+1) = -El(ss+1+m-deltaV1,ii) + Eu(ss+1+m,ii);
                % for second group position of line in cm-1
                %                 DE2(ii,ss+1) = -El(ss+1+n,ii) + Eu(ss+1+n+deltaV2,ii);
            end
            
            n = m + 1;
            t = n + 1;
             for sss = 0:1:peakno+5
                % for first group position of line in cm-1
                DE1s(ii,sss+1) = -El(sss+1+m-deltaV1,ii) + Eu(sss+1+m,ii);
                % for second group position of line in cm-1
                DE2s(ii,sss+1) = -El(sss+1+n-deltaV2,ii) + Eu(sss+1+n,ii);
                
                DE3s(ii,sss+1) = -El(sss+1+t-deltaV3,ii) + Eu(sss+1+t,ii);
                DE4s(ii,sss+1) = -El(sss+1+t-deltaV4,ii) + Eu(sss+1+t,ii);
                %difference between line positin of two groups
%                 DE12(ii,ss+1) = DE2(ii,ss+1)-DE1(ii,ss+1);
            end
            
            
            %difference of transition line position from real data compared to calculated values
            DDE10(ii,:)  = Y1-DE1(ii,:);
            DDE1(ii,:)  = DDE10(ii,:)./Y1;
            %  Least Squares fitting error
            diff22(ii) = sum(DDE1(ii,:).*DDE1(ii,:),2);
            diff33(ii) = sum(DDE10(ii,:).*DDE10(ii,:),2);
            RMSD(ii) = sqrt(diff33(ii)/(peakno));% rmsd error
            %=====================================================================================
            Total(ii,:) = [m,deltaV1,Weuf(ii),WeXeuf(ii),Vef(ii),diff22(ii),ii,i,RMSD(ii)];
            ii=ii+1;
        end
    end
end
Total(all(Total==a,2),:)=[];
Total(all(Total(:,3)<0,2),:)=[];
Total(all(Total(:,4)<0,2),:)=[];
Total(all(Total(:,4)>10,2),:)=[];
Total(all(Total(:,5)<0,2),:)=[];
% stdTotal = std(Total,0,1);

DE1t = transpose(DE1s);
DE2t = transpose(DE2s);
DE3t = transpose(DE3s);
DE4t = transpose(DE4s);

%==================================================================================================================================================
[R2,C2] = ind2sub(size(Total), find(Total(:,9)==min((Total(:,9))))); z2=[R2,C2]; r2=z2(1,1);
fprintf('RMSD minimum is %e, the first index is (Vu = Vl+deltaV)(m deltaV1) %d %d\n', min(Total(:,9)), Total(r2,1:2));
fprintf('values: Weu %4.3f WeXeu %2.6f Ve %6.1f; Wel %4.3f WeXel %2.6f Vel %6.1f ii %d i %d\n',Total(r2,3:5),Wel(Total(r2,8)), WeXel(Total(r2,8)),Vel(Total(r2,8)),Total(r2,7),Total(r2,8));

fprintf('=================== FINAL RESULT =====================\n');
[R2,C2] = ind2sub(size(Total), find(Total(:,6)==min((Total(:,6))))); z2=[R2,C2]; r2=z2(1,1);
fprintf('diff5(total) minimum is %e, the first index is (Vu = Vl+deltaV)(m deltaV1) %d %d\n', min(Total(:,6)), Total(r2,1:2));
fprintf('values: Weu %4.3f WeXeu %2.6f Ve %6.1f; Wel %4.3f WeXel %2.6f Vel %6.1f ii %d i %d\n',Total(r2,3:5),Wel(Total(r2,8)), WeXel(Total(r2,8)),Vel(Total(r2,8)),Total(r2,7),Total(r2,8));
fprintf('======================================================\n');
%===================================================================================================================================================
% toc
