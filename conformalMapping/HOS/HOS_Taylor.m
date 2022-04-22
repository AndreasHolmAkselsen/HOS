function Y_t = HOS_Taylor(t,Y,param)
% Method supports vectorized row input.
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
FFTeta = fft(eta);
FFTvphiS = fft(vphiS);
kx = [0:ceil(N/2)-1, -floor(N/2):-1]';
k = abs(kx);

zz = param.map.xi+1i*eta*param.dim.L;
h = imag(param.map.fz(zz))/param.dim.L;
JInv= abs(param.map.dfz(zz)).^(-2);

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
vphiS_t = - h + JInv.* wNl.*( -.5*vphiS_xi.^2  + .5*(1+eta_xi.^2).*w.^2 );

% Unpad, lowpass filter and dampen:
M = ceil(N/2)-1;
Md = param.kd__kmax*M;
mu = param.rDamping*M*((k-Md)/(M-Md)).^2.*(k>Md);     
kFilter = k<=param.iModeCut;  
Y_t = real([ ifft(kFilter.*fftPad(fft(vphiS_t),N)-mu.*FFTvphiS)
             ifft(kFilter.*fftPad(fft(eta_t  ),N)-mu.*FFTeta  ) ]);        
end

