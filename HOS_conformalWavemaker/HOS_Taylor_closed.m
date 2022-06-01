function Y_t = HOS_Taylor_closed(t,Y,param)

global timeReached
[vphiS,eta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/param.t_end);
end

N = size(eta,1);
dxi = param.map.xi(2)-param.map.xi(1);
dk = pi/(dxi*(N-1));
kx = (0:N-1)'*dk;
k = kx;

% w is the vertical velocity in the zeta-plane; \varphi_\sigma
[w_lin,w_nl] = phiComponentsHOS_closed(vphiS,eta,k,param);
w = w_lin+w_nl;
FFTeta = cosfft(eta);
FFTvphiS = cosfft(vphiS);

zzS     = param.map.xi+1i*eta;
h       = param.map.fy(zzS,t);
JInv    = param.map.fJInv(zzS,t);
ft__fzz = param.map.ft__fzz(zzS,t);


% zzS     = param.map.xi+1i*eta;
% h       = imag(zzS);%param.map.fy(zzS,t);
% JInv    = 1;%param.map.fJInv(zzS,t);
% ft__fzz = param.map.ft__fzz(zzS,t);


if param.DO_PADDING
    error('padding not yet supported in closed domain simulation. Check whether it is straight forward.')
    Nd = N*(4+1)/2;
    w = ifft(fftPad(fft(w),Nd));
    h = ifft(fftPad(fft(h),Nd));
    JInv = ifft(fftPad(fft(JInv),Nd));
    ft__fzz = ifft(fftPad(fft(ft__fzz),Nd));
else
    Nd = N;
end

eta_xi  =  isinfft(fftPad(-kx.*FFTeta,Nd));
vphiS_xi =  isinfft(fftPad(-kx.*FFTvphiS,Nd));

eta_t  =  JInv.*(w.*(1+eta_xi.^2) - vphiS_xi.*eta_xi)  + eta_xi.*real(ft__fzz)-imag(ft__fzz);
vphiS_t = vphiS_xi.*real(ft__fzz) - .5*JInv.*(vphiS_xi.^2-(1+eta_xi.^2).*w.^2 ) - h;

% Unpad, lowpass filter and dampen:
M = N;%ceil(N/2)-1;
Md = param.kd__kmax*M;
mu = param.rDampingDim*M*(((0:N-1)'-Md)/(M-Md)).^2.*((0:N-1)'>Md);    
% kMax = (N-1)*dk;
% kd = param.kd__kmax*kMax*dk;
% mu = param.rDamping*kMax*((k-kd)/(kMax-kd)).^2.*(k>kd);     
kFilter = k<=param.iModeCut*dk;  
Y_t = real([ icosfft(kFilter.*fftPad(cosfft(vphiS_t),N)-mu.*FFTvphiS)
             icosfft(kFilter.*fftPad(cosfft(eta_t  ),N)-mu.*FFTeta  ) ]);        
end

