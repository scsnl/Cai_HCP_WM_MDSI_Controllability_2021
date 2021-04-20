function [Phi, Bv] = get_model_hrf(L,TR,M)


%%%%%%%%%%%%Kernel for HRF Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.5;
[hrf] = spm_hrf(dt);
%d = ceil(length(hrf)/(L));
d = round(TR/dt);
h = hrf(1:d:end)';
h = h(1:L);
h = h./(max(abs(h)));
dh = diff(h);
dh = [0 dh];
h = [h zeros(1,L-length(h))];
dh = [ dh zeros(1,L-length(dh))];
%Bv = spm_orth([h',dh'],'norm');
Bv = [h' dh'];
% for p = 1:size(Bv,2)
%     Bv(:,p) = Bv(:,p) - mean(Bv(:,p));
% end
%Bv = [Bv ones(length(h),1)];
ix = 1;
Phi = [];
for m = 1:M
    Phi1 = zeros(1,M*L);
    Phi2 = zeros(1,M*L);
 %   Phi3 = zeros(1,M*L);
    Phi1(ix:M:M*L) = Bv(:,1)';
    Phi2(ix:M:M*L) = Bv(:,2)';
  %  Phi3(ix:M:M*L) = Bv(:,3)';
    Phi = [Phi;Phi1;Phi2];
    ix = ix + 1;
end