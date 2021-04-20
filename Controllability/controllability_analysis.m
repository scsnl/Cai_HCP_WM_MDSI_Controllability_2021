function [output] = controllability_analysis(input,approach,ofname)

%%%
% function [output] = controllability_analysis(input)
%
% description:

%
% example:
%
% dependencies:
%
% inputs:
%       - input - theta from MDS 
%           - group x subject x condition x roi x roi
%       - approach - 'eigen' default, 'sparsity', or 'connected_components'
%         if using connected components, please be sure to input
%         Theta_normal from MDS output, rather than Theta.
%
% outputs:
%       - output
%           group x subject x condition x roi
% params:
%
% Carlo de los Angeles, 19-12-02

% Editable parameters:
% if doing sparsity based thresholding: 
percent_sparsity = 10;

% Do Not Change Below: 
    if exist(input,'file')
        groups=whos('-file',input);
        data=load(input);
    elseif isstruct(input)
        groups=whos('-struct',output);
        data=input;
    else
        error('Error: unexpected input to function controllability_analysis');
    end

    percent_sparsity = 10;

    for grp_idx=1:length(groups)
        num_subj=size(data.(groups(grp_idx).name),1);
        num_cond=size(data.(groups(grp_idx).name),2);
        num_roi=size(data.(groups(grp_idx).name),3);
        
        for cond_idx =1:num_cond
            
            if strcmp(approach,'sparsity')
                subjectwise_thresholds = [];
                for subj_idx = 1:num_subj    % Subject No. in T test
                    se_grp=squeeze(data.(groups(grp_idx).name)(subj_idx,cond_idx,:,:));
                    zth = getThresholdSparsity(se_grp, percent_sparsity);
                    subjectwise_thresholds = [subjectwise_thresholds zth];
                    c_graph = se_grp .* (abs(se_grp) >= zth); % weighted
                    for m = 1:num_roi
                        [subj_controllability(m),ranks_cond(m)]   = average_controllability_nodewise(c_graph,m);
                    end
                    output.(groups(grp_idx).name)(subj_idx,cond_idx,:)=squeeze(subj_controllability);
                end               
                
                
            elseif strcmp(approach,'eigen')
                subjectwise_thresholds = [];
                for subj_idx = 1:num_subj    % Subject No. in T test
                    se_grp=squeeze(data.(groups(grp_idx).name)(subj_idx,cond_idx,:,:));
                    bestThresholds = getThresholdEigenAnalysis(se_grp);
                    zth = min(bestThresholds);
                    subjectwise_thresholds = [subjectwise_thresholds zth];
                    c_graph = se_grp .* (abs(se_grp) >= zth); % weighted
                    for m = 1:num_roi
                        [subj_controllability(m),ranks_cond(m)]   = average_controllability_nodewise(c_graph,m);
                    end
                    output.(groups(grp_idx).name)(subj_idx,cond_idx,:)=squeeze(subj_controllability);
                end
            elseif strcmp(approach,'connected_components')
                for subj_idx = 1:num_subj
                    se_grp=squeeze(data.(groups(grp_idx).name)(subj_idx,cond_idx,:,:));
                    [G,zth] = get_threshold_connectedComponents(se_grp);
                    c_graph = se_grp .* (abs(se_grp) >= zth);
                    for m = 1:num_roi
                        [subj_controllability(m),ranks_cond(m)]   = average_controllability_nodewise_cc(c_graph,m);
                    end
                    output.(groups(grp_idx).name)(subj_idx,cond_idx,:)=squeeze(subj_controllability);            

                end
            end   
        end
    end
    
    %output_filename=sprintf('mds_controllability_output_%s.mat',approach);
    save(ofname,'-struct','output');

end

function threshold = getThresholdSparsity(network,percent_sparse)

    network = network - diag(diag(network));
    vals = network(:);
    vals = vals(vals ~= 0);
    vals = vals(~isnan(vals));
    sort_vals = sort(abs(vals));
    ix = round(length(vals) * percent_sparse/100);
    try
        threshold = sort_vals(ix);
    catch
        threshold = nan; 
    end
    
end

function bestThresholds = getThresholdEigenAnalysis(network)

    network = network - diag(diag(network));
    vals = abs(network(:));
    vals = vals(vals ~= 0);
    delta = (max(vals) - min(vals))/1000.0;
    thresholds = min(vals):delta:max(vals);
    bestThresholds = [];

    for node = 1:size(network,1)
        th = getThresholdEigenAnalysisNode(network,node,thresholds);
        bestThresholds = [bestThresholds, th];
    end

end

function bestThreshold = getThresholdEigenAnalysisNode(network,node,thresholds)
    flag_uncontrollable = 1;
    bestThreshold = 0.0;
    eps = 10^-20;
    for ii = 1:length(thresholds)
        th = thresholds(ii);
        if flag_uncontrollable
            B = zeros(size(network));
            A = network .* (abs(network) >= th);
            B(node,node) = mean(mean(abs(A)));
            [V,D] = eig(A');
            for j = 1:size(V,2)
                temp = B'*V(:,j);
                if temp'*temp < eps
                    flag_uncontrollable = 0;
                    break
                else
                    bestThreshold = th;
                end
            end
        end
    end
end

function [Gbest, th_best] = get_threshold_connectedComponents(network)
    th_best = 0.0;
    for th = 0.0:0.1:3.0
        thresholded_network = network .* abs(network >= th);
        G = digraph(thresholded_network);
        bins = conncomp(G);
        if length(unique(bins)) == 1
            th_best = th;
            Gbest = G;
        end
    end
    if ~exist('Gbest','var')
        Gbest=nan;
    end
end

function [ controllability,rank_node ] = average_controllability_nodewise( A,node )

    len=length(A);
    A(1:len+1:len*len)=0;

    Ascaled = A./(1+svds(A,1));  %make it Schur stable

    %define C and D matrices
    C = ones(len,1); %read all outputs
    D = zeros(len,1)'; %set feedforward matrix to 0
    ind = find(A ~= 0);
    B = zeros(len,len);

    B(node,node) = mean(mean(abs(A)));

    sys = ss(Ascaled,B,C',D,0.01); %define the discrete time system
    Wc = gram(sys,'c');
    
    %compute average controllability
    controllability = trace(Wc);
    rank_node = rank(Wc * 10^6);
end

function [ controllability,rank_node ] = average_controllability_nodewise_cc( A,node )

len=length(A);
A(1:len+1:len*len)=0;

Ascaled = A./(1+svds(A,1));  %make it Schur stable

C = ones(len,1); %read all outputs
D = zeros(len,1)'; %set feedforward matrix to 0
ind = find(A ~= 0);
B = zeros(len,len);

B(node,node) = mean(mean(abs(A)));

sys = ss(Ascaled,B,C',D,0.01); %define the discrete time system
Wc = gram(sys,'c');
%compute average controllability
controllability = trace(Wc);
rank_node = rank(Wc * 10^6);
end
