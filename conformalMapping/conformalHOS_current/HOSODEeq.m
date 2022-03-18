function [phiS_t,eta_t] = HOSODEeq(t,phiS,eta)
    global nonLinRamp M surfaceMethod x k_cut timeReached t_end H dW
    
    if t-timeReached > 1
        timeReached = floor(t);
        fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/t_end);
    end
    g = 9.81;
    
    switch surfaceMethod
        case 'Taylor'
            wNl = nonLinRamp(t);
            wCurr=1;
            [w_lin,w_nl,phiS_x,eta_x] = phiComponentsHOS(phiS,eta,H,M);
            w = w_lin+w_nl;
            
            dWS = dW(x + 1i*eta);
            Phi_x =  real(dWS);
            Phi_z = -imag(dWS);
            eta_t  =   w_lin + wNl.*(  w_nl          + eta_x.^2.*w - phiS_x.*eta_x ) - wCurr.*(Phi_x.*eta_x-Phi_z);
            phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*w.^2 )...
                - wCurr.*( Phi_x.*phiS_x + .5*Phi_x.^2 + .5*Phi_z.^2 );
%             eta_t  =   w_lin + wNl.*(  w_nl          + eta_x.^2.*w - phiS_x.*eta_x );
%             phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*w.^2 );            

        case 'decayingConformal'  
            nx = numel(x);
            kx = getKx(x); % x->xi
            k = abs(kx);
            if ~isfinite(H), H = realmax; end
            Lsin = -2./(exp(2*kx.*H)-1); Lsin(1) = 1;
            Lcos = 2./(exp(2*kx.*H)+1);
            FFTeta = fft(eta).*(k<k_cut);          
            df =  1 - ifft( kx.*FFTeta.*Lsin);
            if any(real(df) < 0)% downward mapping -> wave breaking
                [eta_t,phiS_t] = deal(nan(size(k)));  return
            end
            JInv = abs(df).^(-2);
            dww = ifft(1i.*kx.*fft(phiS).*(k<k_cut).*Lcos);

            if strcmp(func2str(dW),'@(zz)0') % re-state for efficiency
                FFTb = fft(  -JInv.*imag(dww) ).*(k<k_cut);
                f0 =  1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(nx*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
                tf =  1i*f0 + 1i*ifft(FFTb.*Lsin.*(k~=0));
                eta_t = imag(tf.*df);
                phiS_t = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - g*eta;
            else
                z = x + 1i*ifft(FFTeta.*Lsin,[],1);%+1i*FFTeta(1)/nx;
                dWS = dW(z);
                FFTb = fft(  -JInv.*(imag(dww)+imag(dWS.*df)) ).*(k<k_cut);
                f0 =  1i*sum(imag(FFTb.*conj(FFTeta)).*kx./(nx*sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
                tf =  1i*f0 + 1i*ifft(FFTb.*Lsin.*(k~=0));
                eta_t = imag(tf.*df);
                phiS_t = real(tf).*real(dww) - .5*JInv.*real(dww.^2) - g*eta - real(dWS./conj(df)).*real(dww) - .5*abs(dWS).^2;
                
%                 phiS_t_ = real(dww.*tf) - .5*abs(dww./df+dWS).^2 - g*eta;
%                 max(abs(phiS_t-phiS_t_))    
            end 
            

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
%     df =  1 + 2*fft(kx.*conj(fft(eta)).*(kx>0&k<k_cut)/nx);
%     if any(real(df) < 0)% downward mapping -> wave breaking
%         [eta_t,phiS_t] = deal(nan(size(k)));  return
%     end
%     JInv = abs(df).^(-2);
%     U = conj(-2i*fft(kx.*conj(fft(phiS)).*(kx>0&k<k_cut))/nx);
%     hb = fft(JInv.*imag(U))/nx;
%     tf =  fft(2i*conj(hb).*(kx>0&k<k_cut)) + 1i*hb(1);
%     eta_t = imag(tf.*df);
%     phiS_t = real(tf).*real(U) - .5*JInv.*real(U.^2) - g*eta; % phi_t = "phiS_t"
%     end
