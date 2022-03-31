function z = fConformal(zeta,eta,H,k_cut)
xi = real(zeta);
nx = size(xi,1);
kx = getKx(xi);
k = abs(kx);
if nargin < 4, k_cut = inf; end

FFTeta = fft(eta,[],1);
if isfinite(H)
%     Lsin = -2./(exp(2*kx.*H)-1-2*(k==0));
%     z = zeta + 1i*ifft(FFTeta.*Lsin.*exp(-kx.*imag(zeta)).*(k<k_cut),[],1);
    LsinExp = -2./(exp(kx.*(2*H+imag(zeta)))-(1+2*(k==0)).*exp(kx.*imag(zeta)));
    z = zeta + 1i*ifft(FFTeta.*LsinExp.*(k<k_cut),[],1);
else
%     z =  zeta + 2i*fft(conj(FFTeta).*exp(kx.*imag(zeta).*(kx>0)).*(k<k_cut&kx>0),[],1)+1i*FFTeta(1)/nx;
    z =  zeta + 2i*ifft(FFTeta.*exp(-kx.*imag(zeta).*(kx<0)).*(k<k_cut&kx<0),[],1)+1i*FFTeta(1)/nx;
end

% nx = numel(xi);
% kx = getKx(xi);
% argH = exp(-kx.*H)./sinh(kx.*H); argH(1) = 0;
% 
% FFTeta = fft(eta,[],1);
% if isfinite(H)
%     z = xi - 1i*ifft(FFTeta.*argH.*(abs(kx)<k_cut),[],1)+1i*FFTeta(1)/nx;
% else
%     z =  xi + 2i*fft(conj(FFTeta).*(abs(kx)<k_cut&kx>0),[],1)+1i*FFTeta(1)/nx;
% end