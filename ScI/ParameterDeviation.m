J = jacobian;
r = residual;
[Q,R] = qr(J,0);
mse = sum(abs(r).^2)/(size(J,1)-size(J,2));
Rinv = inv(R);
Sigma = Rinv*Rinv'*mse;
se = sqrt(diag(Sigma));