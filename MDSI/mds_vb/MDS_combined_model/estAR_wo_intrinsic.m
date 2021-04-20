function [B,d,Q] = estAR_wo_intrinsic(Y,u1,Vm,popt)

[N,M] = size(Y);
J = size(Vm,2);
% Get regressors
Xo = [];
Xo = getRegressors(Y,popt);   %Regressors from data Y  
u1 = u1(:);
X = [];    %   %Include External stimuli as a regressor
U2 = [];
for j = 1:J
    u2 = Vm(:,j);
    U2 = u2(popt+1:end) * ones(1,M);
    X1 = Xo.*U2;
    X = [X X1];
end
X = [X u1(popt+1:end)];
Q = zeros(M,M);
B = zeros(M,M,J);
d = zeros(M,1);
for m = 1:M
    y = Y(2:end,m);
    w = X\y;
    ix = 1;
    for j = 1:J
        B(m,1:M,j) = w(ix:ix+M-1);
        ix = ix + M;
    end
    d(m) = w(end);
    Q(m,m) = ((y-X*w)'*(y-X*w))/(N-length(w));
end


function X = getRegressors(Y,p)

[N,M] = size(Y);
X = [];
L = N-p;
for lag = p:-1:1
    for m = 1:M
        x = Y(lag:lag+L-1,m);
        X = [X x];
    end
end