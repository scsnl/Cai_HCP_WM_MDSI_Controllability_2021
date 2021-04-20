function [Phi, Bv] = get_model_hrf_rev(L,TR,M,dt)


%%%%%%%%%%%%Kernel for HRF Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dt = 0.5;
[hrf] = spm_hrf(dt);
hrf = hrf./(max(abs(hrf)));
d_hrf = diff(hrf);
d_hrf = [0;d_hrf];
%d = ceil(length(hrf)/(L));
d = round(TR/dt);
h = hrf(1:d:end)';
h = h(1:L);
%h = h./(max(abs(h)));
%dh = diff(h);
%dh = [0 dh];
dh = d_hrf(1:d:end)';
dh = dh(1:L);
h = [h zeros(1,L-length(h))];
h = h./max(abs(h));
dh = [ dh zeros(1,L-length(dh))];
dh = dh./max(abs(dh));
Bv = [h' dh'];
ix = 1;
Phi = [];
for m = 1:M
    Phi1 = zeros(1,M*L);
    Phi2 = zeros(1,M*L);
    Phi1(ix:M:M*L) = Bv(:,1)';
    Phi2(ix:M:M*L) = Bv(:,2)';
    Phi = [Phi;Phi1;Phi2];
    ix = ix + 1;
end