function Y_t = HOS_Taylor_open(t,Y,param)
% Normalization: t -> t*(L/g)^1/2, (eta,x,y,H) -> (eta,x,y,H)*L, phi -> phi*(L^3*g)^1/2, (p/rho) -> (p/rho)*L*g, k -> k/L
% g = 9.81; L is chosen as domain length/(2*pi) (such that k_j = j)
global timeReached
[vphiS,eta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/param.t_end);
end

wNl = param.nonLinRamp(t);


N = size(eta,1);
dxi = param.map.xi(2)-param.map.xi(1);
dk = 2*pi/(dxi*N);
kx = [0:ceil(N/2)-1, -floor(N/2):-1]'*dk;
k = abs(kx);

% w is the vertical velocity in the zeta-plane; \varphi_\sigma
[w_lin,w_nl] = phiComponentsHOS_open(vphiS,eta,k,param);

FFTeta = fft(eta);
FFTvphiS = fft(vphiS);

zzS = param.map.xi+1i*eta;
h = param.map.fy(zzS);
JInv= param.map.fJInv(zzS);


if param.DO_PADDING
    Nd = N*(4+1)/2;
    w_lin = ifft(fftPad(fft(w_lin),Nd));
    w_nl = ifft(fftPad(fft(w_nl),Nd));
    h = ifft(fftPad(fft(h),Nd));
    JInv = ifft(fftPad(fft(JInv),Nd));
else
    Nd = N;
end

eta_xi  =  ifft(fftPad(1i*kx.*FFTeta,Nd));
vphiS_xi =  ifft(fftPad(1i*kx.*FFTvphiS,Nd));
w = w_lin+w_nl;

eta_t  =  JInv.*( w_lin + wNl.*(  w_nl  + eta_xi.^2.*w - vphiS_xi.*eta_xi ) );
vphiS_t = - h*param.g + JInv.* wNl.*( -.5*vphiS_xi.^2  + .5*(1+eta_xi.^2).*w.^2 );

% Unpad, lowpass filter and dampen:
M = ceil(N/2)-1;
Md = param.kd__kmax*M;
mu = param.rDampingDim*M*(((0:N-1)'-Md)/(M-Md)).^2.*((0:N-1)'>Md);    
% kMax = (N-1)*dk;
% kd = param.kd__kmax*kMax*dk;
% mu = param.rDamping*kMax*((k-kd)/(kMax-kd)).^2.*(k>kd);     
kFilter = k<=param.iModeCut*dk;  
Y_t = real([ ifft(kFilter.*fftPad(fft(vphiS_t),N)-mu.*FFTvphiS)
             ifft(kFilter.*fftPad(fft(eta_t  ),N)-mu.*FFTeta  ) ]);          
end

