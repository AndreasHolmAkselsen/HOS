function z = fConformal(zeta,eta,H,k_cut)
xi = real(zeta);
nx = numel(xi);
kx = getKx(xi);
k = abs(kx);


FFTeta = fft(eta,[],1);
if isfinite(H)
    Lsin = -2./(exp((2*kx.*H))-1-2*(k==0));% argH(1) = 0;
    z = zeta + 1i*ifft(FFTeta.*Lsin.*exp(-kx.*imag(zeta)).*(k<k_cut),[],1);%+1i*FFTeta(1)/nx;
else
    z =  zeta + 2i*fft(conj(FFTeta/nx).*exp(kx.*imag(zeta)).*(k<k_cut&kx>0),[],1)+1i*FFTeta(1)/nx;
end

% nx = numel(xi);
% kx = getKx(xi);
% argH = exp(-kx.*H)./sinh(kx.*H); argH(1) = 0;
% 
% FFTeta = fft(eta,[],1);
% if isfinite(H)
%     z = xi - 1i*ifft(FFTeta.*argH.*(abs(kx)<k_cut),[],1)+1i*FFTeta(1)/nx;
% else
%     z =  xi + 2i*fft(conj(FFTeta/nx).*(abs(kx)<k_cut&kx>0),[],1)+1i*FFTeta(1)/nx;
% end