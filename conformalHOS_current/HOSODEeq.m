function [phiS_t,eta_t] = HOSODEeq(t,phiS,eta)
global nonLinRamp surfaceMethod x k_cut timeReached t_end H dW kx DO_PADDING

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/t_end);
end
g = 9.81;

switch surfaceMethod
    case 'Taylor'
        wNl = nonLinRamp(t);
        wCurr = wNl;
        [w_lin,w_nl] = phiComponentsHOS(phiS,eta);
        
        N = size(eta,1);
        FFTeta = fft(eta);
        dWS = dW(x + 1i*eta);
        if DO_PADDING
            Nd = N*(4+1)/2;
            w_lin = ifft(fftPad(fft(w_lin),Nd));
            w_nl = ifft(fftPad(fft(w_nl),Nd));
            eta = ifft(fftPad(FFTeta,Nd));
            dWS = ifft(fftPad(fft(dWS),Nd));
            kFilter = abs(getKx(x,Nd))<k_cut; 
        else
            Nd = N;
            kFilter = abs(kx)<k_cut; 
        end
        eta_x  =  ifft(fftPad(1i*kx.*FFTeta,Nd));
        phiS_x =  ifft(fftPad(1i*kx.*fft(phiS),Nd));
        w = w_lin+w_nl;
        
        Phi_x =  real(dWS);
        Phi_z = -imag(dWS);
        eta_t  =   w_lin + wNl.*(  w_nl          + eta_x.^2.*w - phiS_x.*eta_x ) - wCurr.*(Phi_x.*eta_x-Phi_z);
        phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*w.^2 )...
            - wCurr.*( Phi_x.*phiS_x + .5*Phi_x.^2 + .5*Phi_z.^2 );
        
        % Lowpass filter and unpad:
%         if mod(Nd,2)==0
%             kFilter = [0:Nd/2-1,Nd/2:-1:1]'*kx(2)<k_cut;
%         else
%             kFilter = [0:(Nd-1)/2, (Nd-1)/2:-1:1]'*kx(2)<k_cut;
%         end

        cla, plot(kx,abs(fft(eta)),kx,abs(kFilter.*fft(eta)),'--');xlim([0,max(kx)]);ylim([0,.2])

        eta_t = real(ifft(fftPad(kFilter.*fft(eta_t),N)));
        phiS_t = real(ifft(fftPad(kFilter.*fft(phiS_t),N)));
        
    case 'Chalikov'
        N = numel(x);
        if ~isfinite(H), H = realmax; end
        Lsin = -2./(exp(2*kx.*H)-1); Lsin(1) = 1;
        Lcos = 2./(exp(2*kx.*H)+1);
        k = abs(kx);  
        FFTeta = fft(eta);
        FFTphiS = fft(phiS);
        z = x + 1i*ifft(FFTeta.*Lsin);
        dWS = dW(z);
        if isempty(k_cut), k_cut = N/2*kx(2); end
        if DO_PADDING
            Nd = (4.5/2)*N;% chalikov
            Nd = floor(Nd/2)*2;
            dWS = ifftPad(fft(dWS),Nd);
            eta = ifft(fftPad(fft(eta),Nd));
        else
            Nd = N;
        end
        df =  1 - ifft(fftPad( kx.*FFTeta.*Lsin,Nd));
        if any(real(df) < 0)% downward mapping -> wave breaking
            [eta_t,phiS_t] = deal(nan(size(k)));  return
        end
        JInv = abs(df).^(-2);
        dww = ifft(fftPad(1i.*kx.*FFTphiS.*Lcos,Nd));
        FFTb = fft(  -JInv.*(imag(dww)+imag(dWS.*df)) );
        assert(mod(N,2)==0);
%         FFTb = [FFTb(1:N/2);0;FFTb(Nd-N/2+2:Nd)]*N/Nd; %FFTb = [FFTb(1:N/2);2*FFTb(Nd-N/2+1);FFTb(Nd-N/2+2:Nd)];      
        FFTb = fftPad(FFTb,N);
        tf0 =  1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(N*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
        tf =  1i*tf0 + 1i*ifft(fftPad(FFTb.*Lsin.*(k~=0),Nd));
        eta_t = imag(tf.*df);
        phiS_t = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - g*eta - real(dWS./conj(df)).*real(dww) - .5*abs(dWS).^2;
            
%         eta_t0 = eta_t;
%         eta_t_ =  real(ifft( fftPad(fft(eta_t),N) ));
%         phiS_t_ = real(ifft(fftPad(fft(phiS_t),N)));
        
        % Unpad and dampen:
        kmax = k(2)*N/2;
        mu = k_cut*sqrt(2*pi*g*kx(2))*N/2 * (2*k/kmax-1).^2.*(k>kmax/2);
        eta_t =  real(ifft( fftPad(fft(eta_t),N) - mu.*FFTeta ));
        phiS_t = real(ifft(fftPad(fft(phiS_t),N) - mu.*FFTphiS));
        
        
                
% % %         cla;plot(x,real(ifftPad(kFilter.*fft(eta_t),N)),(0:Nd-1)*N/Nd*x(2),eta_t,'--')
%         subplot(211), plot(kx,abs(fft(eta_t_)),kx,abs(fft(eta_t)),'--');xlim([0,max(kx)]);ylim([0,2])
% %         subplot(212), plot(kx,abs(FFTeta),'-');xlim([0,max(kx)]);ylim([0,.2])
%         subplot(212), plot(getKx(x,Nd),abs(fft(eta_t0)),'-');xlim([0,inf]);ylim([0,2])
% % %         kxPad = getKx(x,Nd);
% % %         plot(kxPad,abs(fft(eta_t)));xlim([0,max(kxPad)])
        
%         cla;plot(x,real(ifftPad(fft(eta),N)),(0:Nd-1)*N/Nd*x(2),eta,'--');axis equal
%         figure,plot( (0:Nd)/N,(2*(0:Nd)/N-1).^2,(0:Nd)/N,(2*(0:Nd)/N-1).^2.*((0:Nd)>N/2),'--'   )
        


%         N = numel(x);
%         kx = getKx(x); % x->xi
%         k = abs(kx);
%         if ~isfinite(H), H = realmax; end
%         Lsin = -2./(exp(2*kx.*H)-1); Lsin(1) = 1;
%         Lcos = 2./(exp(2*kx.*H)+1);
%         FFTeta = fft(eta).*(k<k_cut);
%         df =  1 - ifft( kx.*FFTeta.*Lsin);
%         if any(real(df) < 0)% downward mapping -> wave breaking
%             [eta_t,phiS_t] = deal(nan(size(k)));  return
%         end
%         JInv = abs(df).^(-2);
%         dww = ifft(1i.*kx.*fft(phiS).*(k<k_cut).*Lcos);
%         
%         if strcmp(func2str(dW),'@(zz)0') % re-state for efficiency
%             FFTb = fft(  -JInv.*imag(dww) ).*(k<k_cut);
%             tf0 =  1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(N*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
%             tf =  1i*tf0 + 1i*ifft(FFTb.*Lsin.*(k~=0));
%             eta_t = imag(tf.*df);
%             phiS_t = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - g*eta;
%         else
%             z = x + 1i*ifft(FFTeta.*Lsin,[],1);%+1i*FFTeta(1)/N;
%             dWS = dW(z);
%             FFTb = fft(  -JInv.*(imag(dww)+imag(dWS.*df)) ).*(k<k_cut);
%             tf0 =  1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(N*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
%             tf =  1i*tf0 + 1i*ifft(FFTb.*Lsin.*(k~=0));
%             eta_t = imag(tf.*df);
%             phiS_t = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - g*eta - real(dWS./conj(df)).*real(dww) - .5*abs(dWS).^2;
%             %                 phiS_t_ = real(dww.*tf) - .5*abs(dww./df+dWS).^2 - g*eta;
%             %                 max(abs(phiS_t-phiS_t_))
%         end
% %         cla;plot(x,real(ifftPad(kFilter.*fft(eta_t),N)),(0:Nd-1)*N/Nd*x(2),eta_t,'--')
% %         cla, plot(kx,abs(fft(eta_t)),k_cut*[1,1],[0,2],'--');xlim([0,max(kx)]);ylim([0,2])

        
        
end

end


% %% validations
% x0_t = mean(real(df.*tf))
% y0_t = mean(imag(df.*tf))
% f = x + 1i*ifft(fft(eta).*Lsin);
% max(abs(  imag(f)-eta  ))./max(abs( eta ))
% w = ifft(            fft(phiS).*Lcos);
% max(abs( real(w) - phiS )) ./max(abs( phiS ))
% max(abs(  imag(tf)-JInv.*imag(U)  ))./max(abs( JInv.*imag(U)  ))

% %% deep water code:
%     df =  1 + 2*fft(kx.*conj(fft(eta)).*(kx>0&k<k_cut)/N);
%     if any(real(df) < 0)% downward mapping -> wave breaking
%         [eta_t,phiS_t] = deal(nan(size(k)));  return
%     end
%     JInv = abs(df).^(-2);
%     U = conj(-2i*fft(kx.*conj(fft(phiS)).*(kx>0&k<k_cut))/N);
%     hb = fft(JInv.*imag(U))/N;
%     tf =  fft(2i*conj(hb).*(kx>0&k<k_cut)) + 1i*hb(1);
%     eta_t = imag(tf.*df);
%     phiS_t = real(tf).*real(U) - .5*JInv.*real(U.^2) - g*eta; % phi_t = "phiS_t"
%     end
