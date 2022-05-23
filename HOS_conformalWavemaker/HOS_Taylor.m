function Y_t = HOS_Taylor(t,Y,param)
% Normalization: t -> t*(L/g)^1/2, (eta,x,y,H) -> (eta,x,y,H)*L, phi -> phi*(L^3*g)^1/2, (p/rho) -> (p/rho)*L*g, k -> k/L, f -> f*L

% Method supports vectorized row input.
global timeReached
[vphiS,eta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));
t_dim = t*param.dim.t;
if t_dim-timeReached > 1
    timeReached = floor(t_dim);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/param.t_end);
end

% w is the vertical velocity in the zeta-plane; \varphi_\sigma
[w_lin,w_nl] = phiComponentsHOS(vphiS,eta,param);
w = w_lin+w_nl;
N = size(eta,1);
FFTeta = fft(eta);
FFTvphiS = fft(vphiS);
kx = [0:ceil(N/2)-1, -floor(N/2):-1]';
k = abs(kx);

zzS_dim = param.map.xi+1i*eta*param.dim.L;
% h = imag(param.map.fz(zzS))/param.dim.L;
% JInv= abs(param.map.dfz(zzS)).^(-2);
h = param.map.fy(zzS_dim,t_dim)/param.dim.L;
JInv = param.map.fJInv(zzS_dim,t_dim);
ft__fz = param.map.ft__fz(zzS_dim,t_dim)/param.dim.U;

if param.DO_PADDING
    Nd = N*(4+1)/2;
    w = ifft(fftPad(fft(w),Nd));
    h = ifft(fftPad(fft(h),Nd));
    JInv = ifft(fftPad(fft(JInv),Nd));
    ft__fz = ifft(fftPad(fft(ft__fz),Nd));
else
    Nd = N;
end

eta_xi  =  ifft(fftPad(1i*kx.*FFTeta,Nd));
vphiS_xi =  ifft(fftPad(1i*kx.*FFTvphiS,Nd));

eta_t  =  JInv.*(w.*(1+eta_xi.^2) - vphiS_xi.*eta_xi)  + eta_xi.*real(ft__fz)-imag(ft__fz);
vphiS_t = vphiS_xi.*real(ft__fz) - .5*JInv.*(vphiS_xi.^2-(1+eta_xi.^2).*w.^2 ) - h;

% Unpad, lowpass filter and dampen:
M = ceil(N/2)-1;
Md = param.kd__kmax*M;
mu = param.rDamping*M*((k-Md)/(M-Md)).^2.*(k>Md);     
kFilter = k<=param.iModeCut;  
Y_t = real([ ifft(kFilter.*fftPad(fft(vphiS_t),N)-mu.*FFTvphiS)
             ifft(kFilter.*fftPad(fft(eta_t  ),N)-mu.*FFTeta  ) ]);        
end

