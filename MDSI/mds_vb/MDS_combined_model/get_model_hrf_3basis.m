function [Phi, Bv] = get_model_hrf_3basis (L,TR,M,dt)

%dt = 1/1000;
d = round(TR/dt);

%-Canonical hrf
[bf p] = spm_hrf(dt);

% add time derivative
%----------------------------------------------------------------------
dp     = 1;
p(6)   = p(6) + dp;
D      = (bf(:,1) - spm_hrf(dt,p))/dp;
bf     = [bf D(:)];
p(6)   = p(6) - dp;

% add dispersion derivative
%------------------------------------------------------------------

dp    = 0.01;
p(3)  = p(3) + dp;
D     = (bf(:,1) - spm_hrf(dt,p))/dp;
bf    = [bf D(:)];

% Orthogonalise and fill in basis function structure
%--------------------------------------------------------------------------
%bf  =  spm_orth(bf);

bf = bf(1:d:end,:);
bf = bf(1:L,:);
bf = [bf;zeros(L-size(bf,1),3)];

Bv = bf;
for k = 1:size(Bv,2)
    %Bv(:,k) = Bv(:,k)./max(abs(Bv(:,k)));
    Bv(:,k) = Bv(:,k)./sqrt((Bv(:,k)'*Bv(:,k)));
end

ix = 1;
Phi = [];
for m = 1:M
    Phi1 = zeros(1,M*L);
    Phi2 = zeros(1,M*L);
    Phi3 = zeros(1,M*L);
    Phi1(ix:M:M*L) = Bv(:,1)';
    Phi2(ix:M:M*L) = Bv(:,2)';
    Phi3(ix:M:M*L) = Bv(:,3)';
    Phi = [Phi;Phi1;Phi2;Phi3];
    ix = ix + 1;
end
end