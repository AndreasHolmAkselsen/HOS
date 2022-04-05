function Y_t = HOSODEeq_mode(t,Y)
% Normalization: t -> t*(L/g)^1/2, (eta,x,y,H) -> (eta,x,y,H)*L, phi -> phi*(L^3*g)^1/2, (p/rho) -> (p/rho)*L*g, k -> k/L
% g = 9.81; L is chosen as domain length/(2*pi) (such that k_j = j)
% Let M be the range of modes; -M<=j<=M. Assume odd number of modes (when including zero)

[FFTphiS,FFTeta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));

global timeReached DO_PADDING chalikov dW

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/chalikov.t_end);
end

[M2,m] = size(FFTphiS);
M = (M2-1)/2;
assert(isscalar(M));
kx = [0:M,-M:-1]';
k = abs(kx);
kd = chalikov.kd__kmax*M;
H = chalikov.H;
if isfinite(H)
    Lsin = -2./(exp(2*kx.*H)-1); Lsin(1) = 1;
    Lcos = 2./(exp(2*kx.*H)+1);
else
    H = realmax;
    Lsin = 2*(kx<0); Lsin(1) = 1; Lcos = Lsin;
end


if DO_PADDING
    N = 4.5*M;% chalikov
%     N=2*M2; % integer*M2 will ensure that spatial points coinside before and after padding.
else
    N = M2;
end
if chalikov.doCurr
    xi = (0:N-1).'/N*2*pi*chalikov.dim.L;
    dWS = dW( xi+1i*ifft(fftPad( FFTeta.*Lsin,N)) ); 
end
% if any(angle(df) < -pi/2), Y_t = nan(2*M2,m);  return; end% downward mapping -> wave breaking
dww = ifft(fftPad(1i.*kx.*FFTphiS.*Lcos,N));
df =  1 - ifft(fftPad( kx.*FFTeta.*Lsin,N));
JInv = abs(df).^(-2);
% if any(JInv>1e+6),     Y_t = nan(2*M2,m);  return; end% downward mapping -> wave breaking
if chalikov.doCurr
    FFTb = fftPad(fft(-JInv.*(imag(dww)+imag(dWS.*df))),M2); % unpad
    tf0 = 1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(M2*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
    tf =  1i*tf0 + 1i*ifft(fftPad(FFTb.*Lsin.*(k~=0),N));
    phiS_t_AA = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - real(dWS./conj(df)).*real(dww) - .5*abs(dWS).^2;
else
    FFTb = fftPad(fft(-JInv.*imag(dww)),M2); % unpad
    tf0 = 1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(M2*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
    tf =  1i*tf0 + 1i*ifft(fftPad(FFTb.*Lsin.*(k~=0),N));
    phiS_t_AA = real(tf).*real(dww) - .5*JInv.*real(dww.^2);
end
eta_t_AA = imag(tf.*df);

% Unpad and dampen:
mu = chalikov.r*M*((k-kd)/(M-kd)).^2.*(k>kd);
FFTeta_t =  fftPad(fft(eta_t_AA ),M2) - mu.*FFTeta ;
FFTphiS_t = fftPad(fft(phiS_t_AA),M2) - FFTeta  - mu.*FFTphiS;



Y_t = [FFTphiS_t;FFTeta_t];





%         if t>.2
% 
%         figure, plot(x,ifft(fftPad(fft(eta_t_AA ),N)))
%         hold on; plot((0:Nd-1)*x(2)*N/Nd,eta_t_AA)
% %             eta_t_AA_unPad =  ifft( fftPad(fft(eta_t_AA),N) );
% %             k = abs(kx);
% %             subplot(211), plot(k,abs(fft(eta_t_AA_unPad)),k,abs(fft(eta_t)),'--',abs(getKx(x,Nd)),abs(fft(eta_t_AA))*N/Nd,':','linewidth',1);xlim([0,inf]);ylim([0,50])
% %             subplot(212), plot(k,abs(FFTeta));xlim([0,inf]);ylim([0,.1])
% %             
%             phiS_t0_unPad =  ifft( fftPad(fft(phiS_t_AA),M2) );
%             if DO_PADDING, kN=abs([0:2*M,-2*M:0])';else, kN=k;end
%             subplot(311), plot(k,abs(fft(phiS_t0_unPad)),k,abs(FFTphiS_t),'--',kN,abs(fft(phiS_t_AA))*M2/N,':',k,mu*1e-4,'-','linewidth',1);xlim([0,inf]);ylim([0,50])
%             subplot(312), plot(k,abs(FFTphiS));xlim([0,inf]);ylim([0,.1])
%             subplot(313), plot(0:M2-1,ifft(FFTphiS),(0:N-1)*M2/N,phiS_t_AA,'--');xlim([0,inf]);%             
% 
%         
%                 subplot(311), plot(k,real(fft(phiS_t0_unPad)),k,real(FFTphiS_t),'--',kN',real(fft(phiS_t_AA))*M2/N,':',k,mu*1e-4,'-','linewidth',1);xlim([0,inf]);ylim([0,50])
%                         
            
            %             
% % %         kxPad = getKx(x,Nd);
% % %         plot(kxPad,abs(fft(eta_t)));xlim([0,max(kxPad)])
%         
% %         cla;plot(x,real(ifftPad(fft(eta),N)),(0:Nd-1)*N/Nd*x(2),eta,'--');axis equal
% %         figure,plot( (0:Nd)/N,(2*(0:Nd)/N-1).^2,(0:Nd)/N,(2*(0:Nd)/N-1).^2.*((0:Nd)>N/2),'--'   )
% %             z = fConformal(x,eta,H,inf);
% %             figure, plot(z)
% 
% %         figure, plot(k,abs(fftPad(fft(eta_t0),N)),'r',k,abs(fft(eta_t)),'b--',abs(getKx(x,Nd)),abs(fft(eta_t0))*N/Nd,':');
%         end
               

