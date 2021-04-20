function [KF,log_lik] = filter_update_mi(A,d,Q,C,R,y,u,x,V_prev,initial)


% Prediction
if initial
    xpred = x + d*u;
      Vpred = V_prev;
else
    xpred = A*x + d*u;
    Vpred = A*V_prev*A' + Q;
end
e = y - C*xpred; % error (innovation)
ss = size(A,1);
S = C*Vpred*C' + R;
Sinv = inv(S);

K = Vpred*C'*Sinv; % Kalman gain matrix
KF.mf = xpred + K*e;
KF.Vf = (eye(ss) - K*C)*Vpred;
KF.Vfj = (eye(ss) - K*C)*A*V_prev;
%Loglikelihood
C1 = C * Vpred * C'+ R;
C1 = 0.5*(C1 + C1');

% if min(diag(C1)) <= 10^-5
%     C1 = C1 + 10^-5 * eye(size(C1,1));
% end
%C1 = C1 + 10^-3 * eye(size(C1,1));
log_lik = log(mvnpdf( y', (C*xpred)', C1));
%log_lik = gaussian_prob(e, zeros(1,length(e)), S, 1);
