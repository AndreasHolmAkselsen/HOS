function Y_t = HOS_Taylor(t,Y)
% Normalization: t -> t*(L/g)^1/2, (eta,x,y,H) -> (eta,x,y,H)*L, phi -> phi*(L^3*g)^1/2, (p/rho) -> (p/rho)*L*g, k -> k/L
% g = 9.81; L is chosen as domain length/(2*pi) (such that k_j = j)
global taylor timeReached t_end dW DO_PADDING dim

[phiS,eta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/(t_end/dim.t));
end

wNl = taylor.nonLinRamp(t);
wCurr = wNl;
[w_lin,w_nl] = phiComponentsHOS(phiS,eta);

N = size(eta,1);
FFTeta = fft(eta);

if DO_PADDING
    Nd = N*(4+1)/2;
    w_lin = ifft(fftPad(fft(w_lin),Nd));
    w_nl = ifft(fftPad(fft(w_nl),Nd));
    eta = ifft(fftPad(FFTeta,Nd));
else
    Nd = N;
end
kx = [0:ceil(Nd/2)-1, -floor(Nd/2):-1]';
kFilter = abs(kx)<taylor.k_cut*dim.L;
x = (0:Nd-1).'/Nd*2*pi*dim.L;
dWS = dW(x + 1i*eta)/dim.U;

eta_x  =  ifft(fftPad(1i*kx.*FFTeta,Nd));
phiS_x =  ifft(fftPad(1i*kx.*fft(phiS),Nd));
w = w_lin+w_nl;

Phi_x =  real(dWS);
Phi_z = -imag(dWS);
eta_t  =   w_lin + wNl.*(  w_nl          + eta_x.^2.*w - phiS_x.*eta_x ) - wCurr.*(Phi_x.*eta_x-Phi_z);
phiS_t = - eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*w.^2 )...
    - wCurr.*( Phi_x.*phiS_x + .5*Phi_x.^2 + .5*Phi_z.^2 );

% Lowpass filter and unpad:
Y_t = real([ ifft(fftPad(kFilter.*fft(phiS_t),N))
             ifft(fftPad(kFilter.*fft(eta_t ),N)) ]);

end

