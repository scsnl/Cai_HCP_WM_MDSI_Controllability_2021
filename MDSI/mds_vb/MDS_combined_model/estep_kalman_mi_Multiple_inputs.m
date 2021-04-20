function [xsmooth,Vsmooth,VVsmooth,sum_log_lik] = estep_kalman_mi_Multiple_inputs(Model,BDS)

M = BDS.M;
state_dim = (BDS.L)*M;
Y = BDS.Y;              %Observations
T = size(Y,2);
%Setup Model parameters for Filtering and Smoothing
%A = Model.A;
A = BDS.A;
if isfield(BDS,'Bm')
    Bm = BDS.Bm;
    Jm = size(Bm,3); % # of modulatory Inputs
end
d = Model.d;
u = BDS.u;    %External Stimuli
if isfield(BDS,'um')
    um = BDS.um;
end
Q = Model.Q;     %State Covariance Matrix
C = Model.C;  %Output Matrix
R = BDS.R;              %Observation Noise Covariance
mo = BDS.mo;               %Initial State mean
Co = BDS.Vo;                  %Initial State Covaraince
F = [eye(M) zeros(M,state_dim-M)];
Psi = [eye(state_dim-M) zeros(state_dim-M,M)];
sum_log_lik = 0;
%%%%%%%% Filtering/Forward Pass %%%%%%%%%%%%%%%%%
for t = 1:T
    if isfield(BDS,'um')    
        Bmt = zeros(M);
        for j = 1:Jm
            Bmt = Bmt + um(t,j).*Bm(:,:,j);
        end
        A1 = A + Bmt;
        At = [A1*F;Psi];
    else
        At = Model.A;
    end 
    if t == 1
        initial = 1;
        xprev = mo;
        Vf_prev = Co;
    else
        initial = 0;
        xprev = KF(t-1).mf;
        Vf_prev = KF(t-1).Vf;
    end
    %fprintf('................ filter update time point: %d \n',t);
    [KF(t),log_lik] = filter_update_mi(At,d,Q,C,R,Y(:,t),u(t),xprev,Vf_prev,initial);
    sum_log_lik = sum_log_lik + log_lik;
end
%%%%%%%% Smoothing/Backward Pass %%%%%%%%%%%%%%%%%
KF(T).ms = KF(T).mf;
KF(T).Vs = KF(T).Vf;
VVsmooth = zeros(state_dim,state_dim,T);
for t = T-1:-1:1
       
    if isfield(BDS,'um')      
        Bmt = zeros(M);
        for j = 1:Jm
            Bmt = Bmt + um(t,j).*Bm(:,:,j);
        end
        A1 = A + Bmt;
        At = [A1*F;Psi];
    else
        At = Model.A;
    end 
    Vs_future = KF(t+1).Vs;
    xfilt = KF(t).mf;
    Vfilt = KF(t).Vf;
    xsmooth_future = KF(t+1).ms;
    Vfj_future = KF(t+1).Vfj;
    Vf_future = KF(t+1).Vf;
    %fprintf('................ smooth  update time point: %d \n',t);
    [xsmooth,Vs,VVsmooth(:,:,t+1)] = smooth_update_mi(xsmooth_future,Vs_future,xfilt,Vfilt,At,d,Q,u(t+1),Vfj_future,Vf_future);
    KF(t).ms = xsmooth; KF(t).Vs = Vs; %KF(t).J = J;
    %VVsmooth(:,:,t+1) = KF(t+1).Vfj + (KF(t+1).Vs - KF(t+1).Vf)*inv(KF(t+1).Vf)*(KF(t+1).Vfj);
end
VVsmooth(:,:,1) = zeros(state_dim,state_dim);
KF(1).Pj = zeros(state_dim);
xsmooth = zeros(state_dim,T);
Vsmooth = zeros(state_dim,state_dim,T);
for t = 1:T
    xsmooth(:,t) = KF(t).ms;
    Vsmooth(:,:,t) = KF(t).Vs;
end
