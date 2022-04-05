function Y_t = HOS_Taylor(t,Y)
% Method supports vectorized row input.
global taylor timeReached t_end DO_PADDING map
g = 9.81;
[vphiS,eta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/(t_end));
end

wNl = taylor.nonLinRamp(t);
% w is the vertical velocity in the zeta-plane; \varphi_\sigma
[w_lin,w_nl] = phiComponentsHOS(vphiS,eta);

N = size(eta,1);
FFTeta = fft(eta);
% kx = [0:ceil(N/2)-1, -floor(N/2):-1]';
kx = getKx(map.xi);
kFilter = abs(kx)<taylor.k_cut;

zz = map.xi+1i*eta;
h = imag(map.fz(zz));
JInv= abs(map.dfz(zz)).^(-2);

if DO_PADDING
    Nd = N*(4+1)/2;
    w_lin = ifft(fftPad(fft(w_lin),Nd));
    w_nl = ifft(fftPad(fft(w_nl),Nd));
    h = ifft(fftPad(fft(h),Nd));
    JInv = ifft(fftPad(fft(JInv),Nd));
else
    Nd = N;
end


eta_xi  =  ifft(fftPad(1i*kx.*FFTeta,Nd));
vphiS_xi =  ifft(fftPad(1i*kx.*fft(vphiS),Nd));
w = w_lin+w_nl;

eta_t  =  JInv.*( w_lin + wNl.*(  w_nl  + eta_xi.^2.*w - vphiS_xi.*eta_xi ) );
vphiS_t = - g*h + JInv.* wNl.*( -.5*vphiS_xi.^2  + .5*(1+eta_xi.^2).*w.^2 );

% Lowpass filter and unpad:
Y_t = real([ ifft(kFilter.*fftPad(fft(vphiS_t),N))
             ifft(kFilter.*fftPad(fft(eta_t ),N)) ]);

end

