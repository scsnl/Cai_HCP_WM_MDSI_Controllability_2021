% This code compares the average controllability measures and their relation
% to in and out degrees across 100 different random networks.
% We choose a system with 11 nodes as it is consistent with the no of nodes
% of the actual system under consideration.


clear all
clc
close all
r1_12=[];
r3_12=[];r4_12=[];r5_12=[];r6_12=[];
r2_12=[];
%generate a random network
for ite=1:100
    ite
    network = rand(11) - rand(11);
    network = network - diag(diag(network));
    vals = abs(network(:));
    vals = vals(vals ~= 0);
    delta = (max(vals) - min(vals))/1000.0;
    thresholds = min(vals):delta:max(vals);
    bestThresholds = [];
    
    
    for node = 1:size(network,1)
        [A,th] = getThresholdEigenAnalysisNode(network,node,thresholds);
        bestThresholds = [bestThresholds, th];
    end
    
    d = abs(eigs(A,1));
    A = A/(d+1);
    %this step ensures stability of the A matrix, which is needed to compute the Grammian
    %set the diagonal elements of A to zero
    A = A -diag(diag(A));
    
    
    len = length(A);
    C = ones(len,len); %read all outputs
    D = zeros(len,len)'; %set feedforward matrix to 0
    B = zeros(len,len);
    Tr = zeros(len,1);
    Tr1 = zeros(len,1);
    Mod = zeros(len,1);
    phi = zeros(len,1);
    Out = zeros(len,1);
    deg = zeros(len,1);
    %the below code will compute average controllability for each node.
    for m = 1:len
        B(m,m) = max(max(A));
        sys = ss(A,B,C',D,-1); %define the discrete time system
        A_mod =A;
        Wc = gram(sys,'c');
        
        rank(ctrb(A,B));
        rank(Wc);
        %compute average controllability
        Tr(m) = trace(Wc);
        A1 =abs(A);
        Inabs(m) = sum(A1(m,:));
        In(m) = sum(A(m,:));
        Outabs(m) = (sum(A1(:,m)));
        Out(m) = (sum(A(:,m)));
        diff(m) = Out(m) - In(m);
        diffabs(m) = abs(abs(Out(m)) - abs(In(m)));
        B(m,m) = 0;
    end
    r1 = corrcoef(Out,Tr);
    r2 = corrcoef(Outabs,Tr);
    r3 = corrcoef(In,Tr);
    r4 = corrcoef(Inabs,Tr);
    r5 = corrcoef(diffabs,Tr);
    r6 = corrcoef(diff,Tr);
    
    r1_12(ite)=r1(1,2);
    r2_12(ite)=r2(1,2);
    r3_12(ite)=r3(1,2);
    r4_12(ite)=r4(1,2);
    r5_12(ite)=r5(1,2);
    r6_12(ite)=r6(1,2);
    
    
end

input_mean=[mean(r1_12),mean(r2_12),mean(r3_12),mean(r4_12),mean(r6_12)];
bar(input_mean);
input_sd=[std(r1_12,0),std(r2_12,0),std(r3_12,0),std(r4_12,0),std(r6_12,0)];
%ylim(clims)
hold on
%errorbar(mean(r1_12),std(r2_12,0),mean(r1_12),std(r2_12,0),'bx')
er=errorbar(1:5,input_mean,input_sd,input_sd);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off

function [A,bestThreshold] = getThresholdEigenAnalysisNode(network,node,thresholds)
flag_uncontrollable = 1;
bestThreshold = 0.0;
eps = 10^-20;
for ii = 1:length(thresholds)
    th = thresholds(ii);
    if flag_uncontrollable
        B = zeros(size(network));
        A = network .* [(abs(network) >= th)];
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

% function [Gbest, th_best] = get_threshold_connectedComponents(network)
%     th_best = 0.0;
%     for th = 0.0:0.1:3.0
%         thresholded_network = network .* abs(network >= th);
%         G = digraph(thresholded_network);
%         bins = conncomp(G);
%         if length(unique(bins)) == 1
%             th_best = th;
%             Gbest = G;
%         end
%     end
%     if ~exist('Gbest','var')
%         Gbest=nan;
%     end
% end
% %%We have generated a thresholded A matrix.
%
% function [ controllability,rank_node ] = average_controllability_nodewise( A,node )
%
%     len=length(A);
%     A(1:len+1:len*len)=0;
%
%     Ascaled = A./(1+svds(A,1));  %make it Schur stable
%
%     %define C and D matrices
%     C = ones(len,1); %read all outputs
%     D = zeros(len,1)'; %set feedforward matrix to 0
%     ind = find(A ~= 0);
%     B = zeros(len,len);
%
%     B(node,node) = mean(mean(abs(A)));
%
%     sys = ss(Ascaled,B,C',D,0.01); %define the discrete time system
%     Wc = gram(sys,'c');
%
%     %compute average controllability
%     controllability = trace(Wc);
%     rank_node = rank(Wc * 10^6);
% end
%
% function [ controllability,rank_node ] = average_controllability_nodewise_cc( A,node )
%
% len=length(A);
% A(1:len+1:len*len)=0;
%
% Ascaled = A./(1+svds(A,1));  %make it Schur stable
%
% C = ones(len,1); %read all outputs
% D = zeros(len,1)'; %set feedforward matrix to 0
% ind = find(A ~= 0);
% B = zeros(len,len);
%
% B(node,node) = mean(mean(abs(A)));
%
% sys = ss(Ascaled,B,C',D,0.01); %define the discrete time system
% Wc = gram(sys,'c');
% %compute average controllability
% controllability = trace(Wc);
% rank_node = rank(Wc * 10^6);
% end
%
%
%
%
%
%
