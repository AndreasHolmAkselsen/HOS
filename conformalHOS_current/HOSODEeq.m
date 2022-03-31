function [phiS_t,eta_t] = HOSODEeq(t,phiS,eta)
global taylor surfaceMethod x timeReached t_end H dW kx DO_PADDING kd__kmax r_damping

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/t_end);
end
g = 9.81;

switch surfaceMethod
    case 'Taylor'
        wNl = taylor.nonLinRamp(t);
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
            kFilter = abs(getKx(x,Nd))<taylor.k_cut; 
        else
            Nd = N;
            kFilter = abs(kx)<taylor.k_cut; 
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
        eta_t = real(ifft(fftPad(kFilter.*fft(eta_t),N)));
        phiS_t = real(ifft(fftPad(kFilter.*fft(phiS_t),N)));
%         cla, plot(kx,abs(fft(eta)),kx,abs(kFilter.*fft(eta)),'--');xlim([0,max(kx)]);ylim([0,.2])

    case 'Chalikov'
        N = numel(x);
        if ~isfinite(H), H = realmax; end
        Lsin = -2./(exp(2*kx.*H)-1); Lsin(1) = 1;
        Lcos = 2./(exp(2*kx.*H)+1);
        k = abs(kx);  
        kmax = k(2)*N/2;
        kd = kd__kmax*kmax;
        FFTeta = fft(eta);%.*(k<kd);
        FFTphiS = fft(phiS);%.*(k<kd);
        z = x + 1i*ifft(FFTeta.*Lsin);
        dWS = dW(z);
        if DO_PADDING
            Nd = 4.5*N;% chalikov
            Nd = floor(Nd/2)*2;
            dWS = ifft(fftPad(fft(dWS),Nd));
            eta = ifft(fftPad(fft(eta),Nd));
        else
            Nd = N;
        end
        df =  1 - ifft(fftPad( kx.*FFTeta.*Lsin,Nd));
%         if any(real(df) < 0)% downward mapping -> wave breaking
        if any(angle(df) < -pi/2)% downward mapping -> wave breaking
            [eta_t,phiS_t] = deal(nan(size(k)));  return
        end
        JInv = abs(df).^(-2);
        dww = ifft(fftPad(1i.*kx.*FFTphiS.*Lcos,Nd));
        FFTb = fft(  -JInv.*(imag(dww)+imag(dWS.*df)) );
        
        
        FFTb = fftPad(FFTb,N);
        tf0 = 1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(N*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
        tf =  1i*tf0 + 1i*ifft(fftPad(FFTb.*Lsin.*(k~=0),Nd));       

        eta_t_AA = imag(tf.*df);
        phiS_t_AA = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - g*eta - real(dWS./conj(df)).*real(dww) - .5*abs(dWS).^2;
        
        % Unpad and dampen:
        mu = r_damping*sqrt(2*pi*g*kx(2))*N/2 * ((k-kd)/(kmax-kd)).^2.*(k>kd);
%         mu = k_cut*sqrt(2*pi*g*kx(2))*N/2 * (2*k/kmax-1).^2.*(k>kmax/2);
        eta_t =  real(ifft(fftPad(fft(eta_t_AA ),N) - mu.*FFTeta ));
        phiS_t = real(ifft(fftPad(fft(phiS_t_AA),N) - mu.*FFTphiS));
        

%         figure, plot(x,ifft(fftPad(fft(eta_t_AA ),N)))
%         hold on; plot((0:Nd-1)*x(2)*N/Nd,eta_t_AA)

%         if t>.18
% %             eta_t_AA_unPad =  ifft( fftPad(fft(eta_t_AA),N) );
% %             k = abs(kx);
% %             subplot(211), plot(k,abs(fft(eta_t_AA_unPad)),k,abs(fft(eta_t)),'--',abs(getKx(x,Nd)),abs(fft(eta_t_AA))*N/Nd,':','linewidth',1);xlim([0,inf]);ylim([0,50])
% %             subplot(212), plot(k,abs(FFTeta));xlim([0,inf]);ylim([0,.1])
% %             
%             phiS_t0_unPad =  ifft( fftPad(fft(phiS_t_AA),N) );
%             k = abs(kx);
%             subplot(311), plot(k,abs(fft(phiS_t0_unPad)),k,abs(fft(phiS_t)),'--',abs(getKx(x,Nd)),abs(fft(phiS_t_AA))*N/Nd,':','linewidth',1);xlim([0,inf]);ylim([0,50])
%             subplot(312), plot(k,abs(fft(phiS)));xlim([0,inf]);ylim([0,.1])
%             subplot(313), plot(x,phiS,x,phiS_t_AA,'--');xlim([0,inf]);
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
% max(abs(  imag(tf) + JInv.*imag(dww)  ))./max(abs( JInv.*imag(dww)  ))

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
