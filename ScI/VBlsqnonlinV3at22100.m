clear all;
clc;

global m n L1 L2 deltaV1
mm = 0;
nn = 0;
kk = 1;

Total = zeros(mm*nn*(2*kk+1),11);
i = 1;

L1 = 10; % L+1 = length of the vector
L2 = 8;

for m = 0:mm;
    for n = 0:nn;
        for deltaV1 = -kk: min(kk,min(m,(n-1)))
            %deltaV1 = -kk:min(kk,min(m,(n-1)))
            x0 = [22066.2520356379,229.748944015261,222.186217676566,0.877004825565658,0.818273937562715];

            %x0 = [ 22065.8 229.325 221.452 0.738 0.657];
            
            
            lb = [ 10000   200     200     0        0  ];
            ub = [ 30000   350     350     2        2  ];
            [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@EPSV3at22100, x0, lb, ub);
            Total(i,:) = [m, n, deltaV1, deltaV1+1, resnorm, x,i];
            
            i = i+1;
        end
    end
end

Total(all(Total==0,2),:)=[];

fprintf('=================== FINAL RESULT ====================\n');
[R2,C2] = ind2sub(size(Total), find(Total(:,5)==min((Total(:,5))))); z2=[R2,C2]; r2=z2(1,1);
fprintf('RMS error minimum is %e, the first index is (Vu = Vl+deltaV)(m n deltaV1 deltaV2) %d %d %d %d\n', min(Total(:,5)), Total(r2,1:4));
fprintf('the corresponding value are: Ve %6.1f Weu Wel %4.3f %4.3f; WeXeu WeXel %2.6f %2.6f  i %d\n',Total(r2,6:10),Total(r2,11));
fprintf('======================================================\n');