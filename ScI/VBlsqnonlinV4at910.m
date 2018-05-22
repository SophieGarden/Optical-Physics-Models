clear all;
clc;

global m n L1 L2 deltaV1
mm = 3; %According to Franck-condon principle, fluorescence should be from V'=0 to
%lower electronic state
nn = 3;
% kk = 1;

Total = zeros(4*nn,11);
i = 1;

L1 = 4; % L+1 = length of the vector
L2 = 6;

x0 = [10952.7866929048,229.779225893534,222.937876534964,0.691037904193324,0.645718834342044];

%x0 = [ 22065.8 229.325 221.452 0.738 0.657];


lb = [ 9000   200     200     0        0  ];
ub = [ 20000   350     350     2        2  ];



for n = 0:nn;
    %the first peak if the strongest peak, hence it must be from v'=0 to
    %v" or from v' to v"=0; boltzmann distribution
    %deltaV1 = -kk:min(kk,min(m,(n-1)))
    
    if n == 0
        m = 0;
        
        for deltaV1 = nn:-1:-nn;
            [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@EPSV4at910, x0, lb, ub);
            Total(i,:) = [m, n, deltaV1, deltaV1+1, resnorm, x,i];
            i = i+1;
        end
        
    else
        m = n-1;
        for deltaV1 = m;
            [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@EPSV4at910, x0, lb, ub);
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