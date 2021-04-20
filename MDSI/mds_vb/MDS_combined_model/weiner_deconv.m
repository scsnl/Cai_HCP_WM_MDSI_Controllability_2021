function Xest = weiner_deconv(Y,h,sig2)

[M,N] = size(Y);
hext = zeros(1,N);
hext(1:length(h)) = h;
Hw = fft(hext);
Xest = zeros(M,N);
for m = 1:M
    y = Y(m,:);
    Yw = fft(y);
    Xest(m,:) = ifft((conj(Hw).*Yw)./(abs(Hw).^2 + sig2(m,m)));
end