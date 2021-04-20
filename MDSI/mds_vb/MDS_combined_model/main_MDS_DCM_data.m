clear all
close all
clc

addpath(genpath(('/mnt/musk2/home/fmri/fmrihome/SPM/spm8')))
addpath('/mnt/musk2/home/sryali/Work/deconv_connectivity/dynamical_systems/MDS_combined_model')
warning('off')

TR = 3;
L = round(32/TR);          % Embedded Dimension
method = 'L1_woi';

load('/mnt/musk2/home/sryali/Work/sparsePartialCorrelation/SPC_DCM/Data/sim1')
result_dir = '/mnt/musk2/home/sryali/Work/Causal_Analysis_DCM/Results_MDS/';
result_fname = 'rest_sim1.mat';
result_fname = strcat(result_dir,result_fname);

st = 1;
S = 1;
for subj = 1:Nsubjects
    subj
    Xo = ts(st:st + Ntimepoints -1,:)';
    st = st + Ntimepoints;
      
    M = Nnodes; % # of Regions
    N = Ntimepoints; % # of Time Points
    %Vm = task_waveforms(1).Vm;
    %J = size(Vm,2);
    J = 1;
    v = zeros(N,1); %External stimulus
    cnt = 1;
    for s = 1:S
        Y = Xo;
        

        Ts = size(Y,2);
        if Ts < N
            Y1(:,1:Ts) = Y;
            Y1(:,Ts+1:N) = repmat(Y(:,Ts),1,length(Ts+1:N));
        else
            Y1 = Y(:,1:N);
        end
        %Vm = task_waveforms(s).Vm;
        Vm = ones(Ts,1);
        
        data(s).Y = Y1;
        data(s).Vm = Vm;
    end
    
    %%%%%%%%%%% Initialization %%%%%%%%%%%%%%
    Vmc = [];
    vc = [];
    Xc = [];
    Yc = [];
    %[Phi, Bv] = get_model_hrf(L,TR,M);
    [Phi, Bv] = get_model_hrf_rev(L,TR,M,1/1000);
    %[Phi,Bv] = get_model_hrf_3basis(L,TR,M);
    for s = 1:S
        
        Vmc = [Vmc;data(s).Vm];
        vc = [vc;v];
        Y = data(s).Y(:,1:N);
        Yc = [Yc Y];
    end
    for m = 1:M
        y = Yc(m,:);
        R(m,m) = var(y);
    end
    R = eye(M);
    % Weiner Deconvolution
    h = Bv(:,1);
    Xest = weiner_deconv(Yc,h,R);
    [Bm,d,Q] = estAR_wo_intrinsic(Xest',vc,Vmc,1);
    A = zeros(M);
    [B,R1] = initialize_B_R_WD(Yc,Xest,Bv',L);
    %%%%%%%%%%%%%%%% Data for KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y = zeros(M,N,S);
    Um = zeros(N,J,S);
    Ue = zeros(N,S);
    for s = 1:S
        Y(:,:,s) = data(s).Y(:,1:N);
        Vm = data(s).Vm;
        Um(:,:,s) = Vm;
        Ue(:,s) = v;
    end
    BDS.Phi = Phi;
    BDS.Phit = Bv';
    BDS.A = A;
    BDS.Bm = Bm;   %Weights due to Modulatory Inputs
    BDS.Q = Q;  %State Covariance
    BDS.d = d;
    BDS.B = B;  %Beta
    BDS.R = R;  %Ouput Covariance
    BDS.L = L; % Embedded Dimension
    BDS.M = M; %# of regions
    BDS.mo = zeros(L*M,1); %Initial state vector mean
    BDS.Vo =  10^-3*eye(L*M);    % Initial Covariance
    BDS.ao = 10^-10; BDS.bo = 10^-10;
    BDS.co = 10^-10; BDS.do = 10^-10;
    BDS.method = method;
    switch BDS.method
        case 'L1'
            if sum(v(:)) ~= 0
                BDS.Alphab = 0*ones(M,(J+1)*M+1); %L1
            else
                BDS.Alphab = 0*ones(M,(J+1)*M); %L1
            end
            BDS.Alphab_op = 10^-0*ones(M,size(Bv,2)); %L1
        case 'L2'
            BDS.Alphab = 0*ones(M,1); %L2
            BDS.Alphab_op = 0*ones(M,1); %L2
        case 'L1_woi'
            BDS.Alphab = 0*ones(M,(J)*M+1); %L1
            BDS.Alphab_op = 0*ones(M,size(Bv,2)); %L1
        case 'L2_woi'
            BDS.Alphab = 0*ones(M,1); %L2
            BDS.Alphab_op = 0*ones(M,1); %L2
    end
    BDS.flag = 0;
    BDS.tol = 10^-4;    %Tolerance for convergence
    BDS.maxIter = 100;   %Max allowable Iterations
    %[BDS,LL] = vb_em_iterations_all_subjs(BDS,Y,Um,Ue);
    [BDS,LL] = vb_em_iterations_combined_par_convergence(BDS,Y,Um,Ue);
    length(LL)
    BDS.Var_A = ones(M);
    [Mapi_normal,Mapm_normal,Mapd_normal] = compute_normalized_stats_Multiple_inputs(BDS.A,BDS.Bm,BDS.Var_A,BDS.Var_Bm,BDS.d,BDS.Var_d);
    
    
    subj_model_parameters(subj).Theta_normal = Mapm_normal;
    subj_model_parameters(subj).Theta = BDS.Theta;
    subj_model_parameters(subj).Cov_mat = BDS.Cov_mat;
    subj_model_parameters(subj).Q = BDS.Q;
    
end
Mask = reshape(abs(net(25,:,:)) > 0,M,M);
save(result_fname,'subj_model_parameters','Mask')

