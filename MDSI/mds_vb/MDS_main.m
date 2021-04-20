function MDS_main(subj_id)

addpath(genpath('/oak/stanford/groups/menon/software/spm8/'));
addpath(genpath('./Library_SR/'));

TR = ;
L = round(32/TR);          % Embedded Dimension
method = 'L2_woi';

project_path = '';
data_path = '';
result_path = '';
if ~exist(result_path)
    mkdir(result_path);
end


% loading time series data
file_name = '';
data1 = load(file_name);
fprintf('load timeseries data %s\n', file_name);

% session number
S = 1;

% loading design matrix
%design_matrix_dir = '';
%design_file_name = '';
%task_wv1 = load([design_matrix_dir design_file_name]);
%fprintf('load design matrix %s\n', design_file_name);
%

%data_subjects{1} = data1;
%data_subjects{2} = data2;
%wv_subjects{1} = task_wv1;
%wv_subjects{2} = task_wv2;

data_subjects{1}.data = data1.data(subj_id);

% setting output file
result_fname = '';

for subj = 1:length(data_subjects{1}.data)
    fprintf('processing subject: %s\n\n', data_subjects{1}.data(subj).subjectID);
    
    M = size(data_subjects{1}.data(subj).timeseries,2);
    
    cnt = 1;
    for s = 1:S
        Vm = data_subjects{1}.data(subj).Vm;
        J = size(Vm,2);
        Y = (data_subjects{1}.data(subj).timeseries)';
        N = size(Y, 2);
        v{s} = zeros(N,1); %External stimulus
        for m = 1:M
            y = Y(m,:);
            y = detrend(y);
            y = y - mean(y);
            y = y./std(y);
            Y(m,:) = y;
        end
        Ts = size(Y,2);
        if Ts < N
            Y1(:,1:Ts) = Y;
            Y1(:,Ts+1:N) = repmat(Y(:,Ts),1,length(Ts+1:N));
        else
            Y1 = Y(:,1:N);
        end
        
        Ts = size(Vm,1);
        if Ts < N
            Vm1(1:Ts,:) = Vm;
            Vm1(Ts+1:N,:) = repmat(Vm(Ts,:),length(Ts+1:N),1);
            % Vm1(Ts+1:N,:) = Vm(Ts,:);
        else
            Vm1 = Vm(1:N,:);
        end
        data(s).Y = Y1;
        data(s).Vm = Vm1;
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
        vc = [vc;v{s}];
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
        Ue(:,s) = v{s};
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
            if sum(v{:}) ~= 0
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
save(result_fname,'subj_model_parameters')

