function Y_t = HOS_Taylor(t,Y,param)
% Normalization: t -> t*(L/g)^1/2, (eta,x,y,H) -> (eta,x,y,H)*L, phi -> phi*(L^3*g)^1/2, (p/rho) -> (p/rho)*L*g, k -> k/L
% g = 9.81; L is chosen as domain length/(2*pi) (such that k_j = j)
global timeReached
[vphiS,eta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));

if t*param.dim.t-timeReached > 1
    timeReached = floor(t*param.dim.t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/param.t_end);
end

wNl = param.nonLinRamp(t*param.dim.t);
% w is the vertical velocity in the zeta-plane; \varphi_\sigma
[w_lin,w_nl] = phiComponentsHOS(vphiS,eta,param);

N = size(eta,1);
FFTeta = cosfft(eta);
FFTvphiS = cosfft(vphiS);

kx = (0:N-1)';
k = kx;

zzS = param.map.xi+1i*eta*param.dim.L;
h = param.map.fy(zzS)/param.dim.L;
JInv= param.map.fJInv(zzS);


if param.DO_PADDING
    error('padding not yet supported in closed domain simulation. Check whether it is straight forward.')
    Nd = N*(4+1)/2;
    w_lin = ifft(fftPad(fft(w_lin),Nd));
    w_nl = ifft(fftPad(fft(w_nl),Nd));
    h = ifft(fftPad(fft(h),Nd));
    JInv = ifft(fftPad(fft(JInv),Nd));
else
    Nd = N;
end

eta_xi  =  isinfft(fftPad(-kx.*FFTeta,Nd));
vphiS_xi =  isinfft(fftPad(-kx.*FFTvphiS,Nd));
w = w_lin+w_nl;

eta_t  =  JInv.*( w_lin + wNl.*(  w_nl  + eta_xi.^2.*w - vphiS_xi.*eta_xi ) );
vphiS_t = - h + JInv.* wNl.*( -.5*vphiS_xi.^2  + .5*(1+eta_xi.^2).*w.^2 );

% Unpad, lowpass filter and dampen:
M = ceil(N/2)-1;
Md = param.kd__kmax*M;
mu = param.rDamping*M*((k-Md)/(M-Md)).^2.*(k>Md);     
kFilter = k<=param.iModeCut;  
Y_t = real([ icosfft(kFilter.*fftPad(cosfft(vphiS_t),N)-mu.*FFTvphiS)
             icosfft(kFilter.*fftPad(cosfft(eta_t  ),N)-mu.*FFTeta  ) ]);        
end

