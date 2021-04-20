function BDS = vbmstep_woi_L1_all_subjs(KS,BDS,Y,um,Ue)

%M-step
%Estimate Parameters Using MLE


%y = BDS.Y;
%u = BDS.u; % External
%um = BDS.um;    %Modulatory
S = size(Y,3); % Number of Subjects
M = size(Y,1); % Number of Regions
T = size(Y,2); % Number of Samples
Jm = size(um,2); %Number of Modulatory Inputs
state_dim = size(KS(1).xsmooth,1);
F = [eye(M) zeros(M,state_dim-M)];
Alphab = BDS.Alphab;
Alphab_op = BDS.Alphab_op;
ao = BDS.ao; bo = BDS.bo;
co = BDS.co; do = BDS.do;
Li = BDS.L;
%Li = 1;
% Compute am, bm, dm and taum

Theta = zeros(M,(Jm+0)*M + 1);
Var_Theta = zeros(M,(Jm+0)*M + 1);
Cov_mat = [];
Q = zeros(M);
T1 = zeros((Jm+0)*M,(Jm+0)*M);
T2 = zeros((Jm+0)*M,1);
T3 = 0;
T4 = zeros(M,(Jm+0)*M);
T5 = zeros(M,1);
T6 = zeros(size(KS(1).Vsmooth(:,:,1)));

LB = 0; %Initialize the lower bound
for s = 1:S
    xsmooth = KS(s).xsmooth;
    Vsmooth = KS(s).Vsmooth;
    VVsmooth = KS(s).VVsmooth;
    for t = Li:T
        Ut = [];
        for j = 1:Jm
            Ut = [Ut;um(t,j,s)*eye(M)];
        end
        T3 = T3 + Ue(t,s)*Ue(t,s);
        Ft = Ut*F;
        Pn_1 = KS(s).Vsmooth(:,:,t-1) + xsmooth(:,t-1)*xsmooth(:,t-1)';
        T1 = T1 + Ft * Pn_1 * Ft';
        T2 = T2 + Ue(t,s)*Ft * xsmooth(:,t-1);
        Pj = VVsmooth(:,:,t) + xsmooth(:,t)*xsmooth(:,t-1)';
        T4 = T4 + Pj(1:M,1:state_dim)*Ft';
        T5 = T5 + xsmooth(1:M,t)*Ue(t,s);
        Pn = Vsmooth(:,:,t) + xsmooth(:,t)*xsmooth(:,t)';
        T6 = T6 + Pn;

        %Compute Entropy of state variables s(t)
        cov_st = Vsmooth(1:M,1:M,t);
        LB = LB + 0.5 * log(det(cov_st)) + 0.5 *M*( 1 + log(2*pi));

    end
end
for s = 1:S
    Vsmooth = KS(s).Vsmooth;
    for t = 1:1
        cov_st = Vsmooth(1:M,1:M,t);
        LB = LB + 0.5 * log(det(cov_st)) + 0.5 *M*( 1 + log(2*pi));
    end
end

for m = 1:M

    t4 = T4(m,:);
    t5 = T5(m);
    t6 = T6(m,m);

    inv_VN = [T1 T2;T2' T3] + diag(Alphab(m,:));
    VN = pinv(inv_VN);
    wbar = VN * [t4';t5];
    Theta(m,:) = wbar';         %Mean of each parameter
    Var_Theta(m,:) = diag(VN)'; %Variance of each parameter
    aN = ao + 0.5*(S*(T-Li + length(wbar)));
    bN = bo + 0.5 * (t6 - wbar' * inv_VN * wbar);
    cN = co + 0.5;
    for n = 1:length(wbar)
        dN = do + 0.5*( wbar(n)^2 *( aN/bN) + VN(n,n));
        Alphab(m,n) = cN/dN;
    end
    Q(m,m) = bN/aN;
    Var_Theta(m,:) = bN/aN .*diag(VN)'; %Variance of each parameter
    %inv_Var_Theta(m,:) = diag(inv_VN).*aN/bN;
    Cov_mat(:,:,m) = bN/aN .* VN;
    
    
     % Compute KL divergence between p(A,C,d,Q,alpha) &  q((A,C,d,Q,alpha)
    D = length(wbar);
    term1 = 0.5*D*(psi(1,aN) - log(bN) + psi(1,cN) - log(dN)-log(2*pi))  ...
        -0.5*(cN/dN)*((aN/bN)*wbar'*wbar + trace(VN)) - log(gamma(ao)) + ...
        ao*log(bo) + (ao-1)*(psi(1,aN)-log(bN)) - bo * (aN/bN);
    term2 = -log(gamma(co)) + co*log(do) + (co-1)*(psi(1,cN)-log(dN)) - do*(cN/dN);
    term3 = 0.5*D*(psi(1,aN)-log(bN)-log(2*pi)-1)-0.5*log(det(VN)) - log(gamma(aN)) ...
        + aN*log(bN) + (aN-1)*(psi(1,aN) - log(bN)) - aN;
    term4 = - log(gamma(cN)) + (cN-1) * psi(1,cN) + log(dN) - cN;
    LB = LB + ( term1 + term2 - term3 - term4);

end
idx = 1;
Bm = zeros(M,M,Jm);
Var_Bm = zeros(M,M,Jm);
for j = 1:Jm
    Bm(1:M,1:M,j) = Theta(1:M,idx:idx+M-1);
    Var_Bm(1:M,1:M,j) = Var_Theta(1:M,idx:idx+M-1);
    idx = idx + M;
end
d = Theta(:,end);
Var_d = Var_Theta(:,end);
%%%%%%%%%% Estimates of mo and Vo
mo = zeros(size(KS(1).xsmooth(:,1)));
for s = 1:S
    mo = mo + KS(s).xsmooth(:,1);
end
mo = mo./S;
Vo = zeros(size(KS(1).Vsmooth(:,:,1)));
for s = 1:S
    xs1 = KS(s).xsmooth(:,1);
    Vo = Vo + KS(s).Vsmooth(:,:,1) + (mo - xs1)*(mo - xs1)';
end
Vo = Vo./S;
% Estimation of B and R %%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_Pn = zeros(state_dim);
sum_yn = zeros(M);
y1 = Y(:,:,1);          %For initialization
xsmooth = KS(1).xsmooth; %For initialization
for m = 1:M
    temp = y1(m,1)*xsmooth(m:M:end,1); %sum of outerproduct of yn and mn
    sum_yn_mn1{m} = zeros(size(temp));
end
for s = 1:S
    xsmooth = KS(s).xsmooth;
    Vsmooth = KS(s).Vsmooth;
    y = Y(:,:,s);
    for t = Li:T
        sum_Pn = sum_Pn + Vsmooth(:,:,t) + xsmooth(:,t)*xsmooth(:,t)';
        for m = 1:M
            sum_yn_mn1{m} = sum_yn_mn1{m} + y(m,t)*xsmooth(m:M:end,t); %sum of outerproduct of yn and mn
        end
        sum_yn = sum_yn + y(:,t) * y(:,t)';
    end
end

Phi1 = BDS.Phit;
R = zeros(M);
Lavg = 0;
for m = 1:M
    sum_Pn1 = sum_Pn(m:M:end,m:M:end);

    inv_VN = Phi1 * sum_Pn1 * Phi1' + diag(Alphab_op(m,:));
    VN = pinv(inv_VN);
    wbar = VN * Phi1 * sum_yn_mn1{m};
    aN = ao + 0.5*(S*(T-Li + length(wbar)));
    bN = bo + 0.5 * (sum_yn(m,m) - wbar' * inv_VN * wbar);
    cN = co + 0.5;

    for n = 1:length(wbar)
        dN = do + 0.5*( wbar(n)^2 *( aN/bN) + VN(n,n));
        Alphab_op(m,n) = cN/dN;
    end
    R(m,m) = bN/aN;
    B1(:,m) = wbar;
    
    
    % Compute KL divergence between p(A,C,d,Q,alpha) &  q((A,C,d,Q,alpha)
    D = length(wbar);
    term1 = 0.5*D*(psi(1,aN) - log(bN) + psi(1,cN) - log(dN)-log(2*pi))  ...
        -0.5*(cN/dN)*((aN/bN)*wbar'*wbar + trace(VN)) - log(gamma(ao)) + ...
        ao*log(bo) + (ao-1)*(psi(1,aN)-log(bN)) - bo * (aN/bN);
    term2 = -log(gamma(co)) + co*log(do) + (co-1)*(psi(1,cN)-log(dN)) - do*(cN/dN);
    term3 = 0.5*D*(psi(1,aN)-log(bN)-log(2*pi)-1)-0.5*log(det(VN)) - log(gamma(aN)) ...
        + aN*log(bN) + (aN-1)*(psi(1,aN) - log(bN)) - aN;
    term4 = - log(gamma(cN)) + (cN-1) * psi(1,cN) + log(dN) - cN;
    
    LB = LB + ( term1 + term2 - term3 - term4);
    
    Lavg = Lavg + (-T1/2)*log(2*pi) - 0.5*T1*(psi(1,aN)-log(bN)) - 0.5*(aN/bN)* ...
        (sum_yn(m) +  trace( (wbar*wbar' + (bN/aN) * VN) * Phi1 * sum_Pn1 * Phi1') - ...
        2 * wbar' * Phi1 * sum_yn_mn1{m});    
    
    
end
P = size(Phi1,1);
B = zeros(M,M*P);
ix = 1;
for m = 1:M
    B(m,ix:ix+P-1) = B1(:,m)';
    ix = ix + P;
end

A = zeros(M);
%d = zeros(M,1);
BDS.A = A;
BDS.Bm = Bm;
BDS.d = d;
BDS.Q = Q;
BDS.B = B;
BDS.R = R;
BDS.Alphab = Alphab;
BDS.Alphab_op = Alphab_op;
BDS.Var_Bm = Var_Bm;
BDS.Var_d = Var_d;
BDS.Cov_mat = Cov_mat;
BDS.Theta = Theta;
BDS.mo = mo;
BDS.Vo = Vo;
BDS.LB = LB;


diag_R = diag(R);
neg_elements = diag_R < 0;
if sum(neg_elements) > 0
    BDS.flag = 1;
end