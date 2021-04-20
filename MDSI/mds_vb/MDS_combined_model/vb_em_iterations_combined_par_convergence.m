function [BDS,LL,KS] = vb_em_iterations_combined_par_convergence(BDS,Y,Um,Ue)

%Setup Parameters for Kalman Smoothing (E-step);
M = size(BDS.A,1);
state_dim = (BDS.L)*M;
F = [eye(M) zeros(M,state_dim-M)];
Psi = [eye(state_dim-M) zeros(state_dim-M,M)];
T = size(Y,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Iterate between E and M steps Until Convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
previous_LB = -inf;
converged = 0;
iter = 1;
LL = zeros(BDS.maxIter,1);
Theta_prev = [BDS.A(:); BDS.Bm(:); BDS.d(:); BDS.Q(:); BDS.R(:)];
while ~converged && (iter <= BDS.maxIter) 

    %iter
    %Setup Model parameters for Filtering and Smoothing
    Model.A = [(BDS.A)*F;Psi];   %State Transition Natrix
    d = zeros(state_dim,1);
    d(1:M) = BDS.d;       %Input Vector
    Model.d = d;
    Q = zeros(state_dim);
    Q(1:M,1:M) = BDS.Q;     %State Covariance Matrix
    Model.Q = Q;
    Model.C = (BDS.B)*(BDS.Phi);  %Output Matrix
    %%%%% E step
    BDS_old = BDS;
    for s = 1:size(Y,3)
        BDS.Y = Y(:,:,s);
        BDS.um = Um(:,:,s);
        BDS.u = Ue(:,s);
%         if iter > 1
%             BDS.mo = KS(s).xsmooth(:,1);
%             BDS.Vo = KS(s).Vsmooth(:,:,1);
%         end
        fprintf('-------E-step kalman filter at session : %d\n', s);
        [KS(s).xsmooth,KS(s).Vsmooth,KS(s).VVsmooth,sub_loglik(s)] = estep_kalman_mi_Multiple_inputs(Model,BDS);
        %BDS.xs = xsmooth(1:M,1:T);
    end
    %%%%%%%%%% M step
    display('------- M-step -----');
    if ~isfield(BDS,'um');
        BDS = vbmstep_ext_inputs(KS,BDS);
    else
        switch BDS.method
            case 'L1'
                %BDS = vbmstep_Multiple_mod_inputs(xsmooth,Vsmooth,VVsmooth,BDS);
                BDS = vbmstep_L1(xsmooth,Vsmooth,VVsmooth,BDS);
            case 'L2'
                BDS = vbmstep_L2(xsmooth,Vsmooth,VVsmooth,BDS);
            case 'L1_woi'
                BDS = vbmstep_woi_L1_all_subjs(KS,BDS,Y,Um,Ue);
            case 'L2_woi'
                BDS = vbmstep_woi_L2_all_subjs(KS,BDS,Y,Um,Ue);
        end
    end
    sub_loglik = sub_loglik(abs(sub_loglik) ~= Inf);
    LB = sum(sub_loglik);
    Theta = [BDS.A(:); BDS.Bm(:); BDS.d(:); BDS.Q(:); BDS.R(:)];
    %Check for convergence
    err = ((Theta(:) - Theta_prev(:))'*(Theta(:) - Theta_prev(:))/(Theta(:)'*Theta(:)));
    fprintf('At iteration: %d, err: %d , LB: %d \n',iter, err, LB);
    %err
    if err <= BDS.tol
        converged = 1;

    else
        Theta_prev = Theta;
        LL(iter)=LB;
        iter = iter + 1;
        %LL = [LL LB];
        
        
    end
    

end













function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
% EM_CONVERGED Has EM converged?
% [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
%
% We have converged if the slope of the log-likelihood function falls below 'threshold',
% i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
% 'threshold' defaults to 1e-4.
%
% This stopping criterion is from Numerical Recipes in C p423
%
% If we are doing MAP estimation (using priors), the likelihood can decrase,
% even though the mode of the posterior is increasing.

if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
    if loglik - previous_loglik < -1e-3 % allow for a little imprecision
        fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
        decrease = 1;
        converged = 0;
        return;
    end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
%avg_loglik = abs(loglik);
if (delta_loglik / avg_loglik) < threshold, converged = 1; end



