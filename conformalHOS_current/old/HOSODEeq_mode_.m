function [FFTphiS_t,FFTeta_t] = HOSODEeq_mode(t,FFTphiS,FFTeta)
global x timeReached t_end H dW kx DO_PADDING kd__kmax r_damping

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/t_end);
end
g = 9.81;


N = numel(x);
if ~isfinite(H), H = realmax; end
Lsin = -2./(exp(2*kx.*H)-1); Lsin(1) = 1;
Lcos = 2./(exp(2*kx.*H)+1);
k = abs(kx);
kmax = k(2)*N/2;
kd = kd__kmax*kmax;
z = x + 1i*ifft(FFTeta.*Lsin);
dWS = dW(z);
if DO_PADDING
    Nd = 4.5*N;% chalikov
    Nd = floor(Nd/2)*2;
    dWS = ifft(fftPad(fft(dWS),Nd));
else
    Nd = N;
end
%         eta = ifft(fftPad(FFTeta,Nd));
df =  1 - ifft(fftPad( kx.*FFTeta.*Lsin,Nd));
%         if any(real(df) < 0)% downward mapping -> wave breaking
if any(angle(df) < -pi/2)% downward mapping -> wave breaking
    [FFTeta_t,FFTphiS_t] = deal(nan(size(k)));  return
end
JInv = abs(df).^(-2);
dww = ifft(fftPad(1i.*kx.*FFTphiS.*Lcos,Nd));
FFTb = fft(  -JInv.*(imag(dww)+imag(dWS.*df)) );


FFTb = fftPad(FFTb,N);
tf0 = 1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(N*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
tf =  1i*tf0 + 1i*ifft(fftPad(FFTb.*Lsin.*(k~=0),Nd));

eta_t_AA = imag(tf.*df);
phiS_t_AA = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - real(dWS./conj(df)).*real(dww) - .5*abs(dWS).^2;

% Unpad and dampen:
mu = r_damping*sqrt(2*pi*g*kx(2))*N/2 * ((k-kd)/(kmax-kd)).^2.*(k>kd);
%         mu = k_cut*sqrt(2*pi*g*kx(2))*N/2 * (2*k/kmax-1).^2.*(k>kmax/2);
FFTeta_t =  fftPad(fft(eta_t_AA ),N) - mu.*FFTeta ;
FFTphiS_t = fftPad(fft(phiS_t_AA),N)  - g*FFTeta  - mu.*FFTphiS;




end

