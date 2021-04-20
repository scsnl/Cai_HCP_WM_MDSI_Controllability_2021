function [B,R,B1] = initialize_B_R_WD(Y,S,Phi,L)

% Formulate Vt
M = size(Y,1);
N = size(Y,2);

P = size(Phi,1);     % # of regressors
for m = 1:M
    Mat = zeros(P);
    Rs = zeros(P,1);
    X = [];
    for n = L+1:N
        sn = S(m,n:-1:n-L+1)';
        xn = Phi * sn;
        Mat = Mat + xn * xn';
        Rs = Rs + xn*Y(m,n);
        X = [X;xn'];
    end
    %b = pinv(Mat)*Rs;
    b = Mat\Rs;
    B1(:,m) = b;
    % Variance
    y = Y(m,L+1:N)';
    R(m,m) = ((y - X*b)'*(y - X*b))/(N-L);
end
B = zeros(M,M*P);
ix = 1;
for m = 1:M
    B(m,ix:ix+P-1) = B1(:,m)';
    ix = ix + P;
end