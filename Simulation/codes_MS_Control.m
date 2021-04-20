clear all
clc
close all
%generate a random network of 'n' nodes. The adjacency matrix can have
%signed weights.
n = 100,;
network = rand(n) - rand(n);
%The next line sets the diagonal elements to zero
network = network - diag(diag(network));

%The random network generated will result in a complete graph, which may
%not be useful for our analysis. We threshold the graph in such a way that
%renders the system controllable by each node.
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
%the previous step ensures stability of the A matrix, which is needed to compute the Grammian


%we next compute the trace of the Grammian Matrix and compare it with
%different graph properties, like the in weights, out weights and the
%difference.
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

%r1-r6 compute the correlation coefficients between the Average
%Controllability and several degree measures.
r1 = corrcoef(Out,Tr);
r2 = corrcoef(Outabs,Tr);
r3 = corrcoef(In,Tr);
r4 = corrcoef(Inabs,Tr);
r5 = corrcoef(diffabs,Tr);
r6 = corrcoef(diff,Tr);

%The first figure here gives a plot between the Out degrees vs the Average
%Controllability where as the second figure is a plot between absolute
%weigted out degree and the Average Controllability. These plots appear in Supplementary Figure 2
figure('Name', 'Average Controllability vs Out Degrees')
subplot(121)
plot(Out,Tr,'o')
xlabel('Weighted Out degree')
ylabel('Average Controllability')
title(r1(1,2))
subplot(122)
plot(Outabs,Tr,'*')
xlabel('Absolute of weighted Out degree')
ylabel('Average Controllability')
title(r2(1,2))
%

%The first figure here gives a plot between the In degrees vs the Average
%Controllability where as the second figure is a plot between absolute
%weigted In degree and the Average Controllability. These plots appear in Supplementary Figure 3
figure('Name', 'Average Controllability vs In Degrees')
subplot(121)
plot(In,Tr,'o')
xlabel('Weighted In degree' )
ylabel('Average Controllability')
title(r3(1,2))
subplot(122)
plot(Inabs,Tr,'*')
xlabel('Absolute of Weighted In degree')
ylabel('Average Controllability')
title(r4(1,2))

% %The first figure gives a plot between the difference in Out and In degrees vs the Average
%Controllability and the second figure plots the Average Controllability vs
%the absolute values on the Out-In degrees. These plots appear in
%Supplementary Figure 4.
%Note: In some randomly generated networks, one may find spurious correlations
%between the average controllability and the absolute value of Out-In
%degrees. The variance of the correlations is very high when one repeats
%the experiment multiple times. This is shown by using the script
%"average_100subjects.m" in Supplementary Figure-5.
figure('Name', 'Average Controllability vs Out-In Degrees')
subplot(122)
plot(diffabs,Tr,'o')
xlabel('Absolute value of Out-In Degree' )
ylabel('Average Controllability')
title(r5(1,2))
subplot(121)
plot(diff,Tr,'*')
xlabel('Out-In degree degree')
ylabel('Average Controllability')
title(r6(1,2))


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

%