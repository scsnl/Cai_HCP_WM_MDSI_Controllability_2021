function [outflow,inflow,net_outflow] = thresholdedNetworkForDegreeAnalysis(subj_model_parameters,pval,condition,threshold)


% The script finds the weighted outflow, weighted inflow, and the net
% weighted outflow from the MDSI estimated causal network for a given task
% condition

%Inputs
%subj_model_parameters - MDS estimated causal matrices
%pval : p-value at which network needs to be thresholded (default: 
%       use pval= 0.05
% condition: Task condition
%Threshold: 1 - Apply threshold at pval
%           0 - No thhreshold

%Output
%outflow :  weighted out flow
%inflow :  weighted in flow
%net_outflow: weighted out flow-weighted in flow


no_subjs = length(subj_model_parameters);
if isempty(pval)
    pval = 0.05;
end

for subj = 1:no_subjs    % Subject No. in T test
    % Raw MDS weights
    se_grp = getWeights(subj_model_parameters(subj), condition);
    % Normal MDS weights
    wts_normal = subj_model_parameters(subj).Theta_normal(:,:,condition);
    zth = abs(norminv(pval));
    binary_graph = abs(wts_normal) >= zth; %weighted
    
    if threshold
        c_graph = se_grp .* binary_graph  ; %weighted
    else
        c_graph = se_grp;
    end
    
    outflow(subj,:) = nansum(c_graph, 1); % using nansum here ignores diagonal entries
    inflow(subj,:)  = nansum(c_graph, 2)'; % flip to make row vector
    net_outflow(subj,:)   = outflow(subj,:) - inflow(subj,:);
end

function Theta = getWeights(subj_model_parameters, condition)
M = size(subj_model_parameters(1).Theta,1); % No of nodes
for subj = 1:length(subj_model_parameters)
    Theta = subj_model_parameters(subj).Theta(:, (condition-1)*M+1:condition*M);
end




