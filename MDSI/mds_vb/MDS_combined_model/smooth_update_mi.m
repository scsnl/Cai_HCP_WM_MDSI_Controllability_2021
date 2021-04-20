function [xsmooth,Vsmooth,VVsmooth_future] = smooth_update_mi(xsmooth_future,Vsmooth_future,xfilt,Vfilt,A,d,Q,u,VVfilt_future,Vfilt_future)

xpred = A*xfilt + d*u;
Vpred = A*Vfilt*A' + Q; % Vpred = Cov[X(t+1) | t]
J = Vfilt * A' * inv(Vpred); % smoother gain matrix
xsmooth = xfilt + J*(xsmooth_future - xpred);
Vsmooth = Vfilt + J*(Vsmooth_future - Vpred)*J';
VVsmooth_future = VVfilt_future + (Vsmooth_future - Vfilt_future)*inv(Vfilt_future)*VVfilt_future;