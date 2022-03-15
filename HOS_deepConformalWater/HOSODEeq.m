function [phiS_t,eta_t] = HOSODEeq(t,phiS,eta)
    global nonLinRamp M surfaceMethod x k_cut timeReached t_end H

    if t-timeReached > 1
        timeReached = floor(t);
        fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/t_end);
    end
    wNl = nonLinRamp(t); 
    g = 9.81;
    
    switch surfaceMethod
        case 'Taylor'
            [w_lin,w_nl,phiS_x,eta_x] = phiComponentsHOS(phiS,eta,H,M);
            W = w_lin+w_nl;
            eta_t  =   w_lin + wNl.*(  w_nl          + eta_x.^2.*W - phiS_x.*eta_x );
            phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*W.^2 );            
        case {'phiS_x','phi_x'}
            [w_lin,~,phiS_x,eta_x] = phiComponentsHOS(phiS,eta,H,1);
            U = phiComponentsConformal(phiS,eta);
            w = imag(U); u = real(U);
            w_nl = w-w_lin;
            switch surfaceMethod
                case 'phiS_x'
                    eta_t  =   w_lin + wNl.*(  w_nl          + eta_x.^2.*w - phiS_x.*eta_x );
                    phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*w.^2 );
                case 'phi_x'
                    eta_t  = w_lin +  wNl.*(w_nl - u.*eta_x);
                    phiS_t = - g*eta + wNl.*(eta_t.*w -.5*(u.^2+w.^2));
            end
        case 'decayingConformal'  
            nx = numel(x);
            kx = getKx(x); % x->xi
            k = abs(kx);
             
            if isfinite(H)
                argH = 2./(exp((2*kx.*H))-1).*(k<k_cut); argH(1) = 0;
                heta = fft(eta)/nx;
                df =  1 + ifft( kx.*heta.*argH)*nx;
                if any(real(df) < 0)% downward mapping -> wave breaking
                    [eta_t,phiS_t] = deal(nan(size(k)));  return
                end
                JInv = abs(df).^(-2);
                U = conj(ifft(-1i.*kx.*fft(phiS).*argH));
                hb = fft(JInv.*imag(U))/nx;
                f0 = hb(1)+1i*sum(imag(hb.*conj(heta)).*kx./(sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
                tf =  1i*f0 -1i*ifft(hb.*argH)*nx;
                
%                 x0_t = mean(real(df.*tf))
% %                 f = x - 1i*ifft(fft(eta).*argH.*(k<k_cut))+1i*mean(eta); 
% %                 figure;plot(x,imag(f),'-',x,eta,'--','linewidth',1.5)
%                 max(abs(  imag(f)-eta  ))./max(abs( eta ))
%                 w = ifft(            -fft(phiS).*argH.*(k<k_cut))+mean(phiS); 
% %                 figure;plot(x,real(w),'-',x,phiS,'--','linewidth',1.5)
%                 max(abs( real(w) - phiS )) ./max(abs( phiS ))
% %                 figure;plot(x,imag(tf),'-',x,JInv.*imag(U),'--','linewidth',1.5)
%                 max(abs(  imag(tf)-JInv.*imag(U)  ))./max(abs( JInv.*imag(U)  ))


            else
                df =  1 + 2*fft(kx.*conj(fft(eta)).*(kx>0&k<k_cut)/nx);
                if any(real(df) < 0)% downward mapping -> wave breaking
                    [eta_t,phiS_t] = deal(nan(size(k)));  return
                end
                JInv = abs(df).^(-2);
                U = conj(-2i*fft(kx.*conj(fft(phiS)).*(kx>0&k<k_cut))/nx);
                hb = fft(JInv.*imag(U))/nx;
                tf =  fft(2i*conj(hb).*(kx>0&k<k_cut)) + 1i*hb(1);
            end
            eta_t = real(df).*imag(tf) + imag(df).*real(tf);
            phiS_t = real(tf).*real(U) - .5*JInv.*real(U.^2) - g*eta; % phi_t = "phiS_t"
            
%             
%             max(abs(df_-df))/max(abs(df_))
%             max(abs(U_-U))/max(abs(U_))
%             max(abs(tf_-tf))/max(abs(tf_))
            
    end

    
%     heta_t = fft(eta_t); hphiS_t = fft(phiS_t);     
%     eta_t  = ifft(heta_t.*(k>0 & k<k_cut_ & abs(heta_t)>1e-6));
%     phiS_t = ifft(hphiS_t.*(k>0 & k<k_cut_ & abs(hphiS_t)>1e-6));     
         
end

     %%%      adding a source term.
%             if max(abs(eta))< etaCap               
%                 phiS_t = phiS_t + sourceAmplification *phiS/max(abs(phiS));
%                 eta_t = eta_t + sourceAmplification *eta_t/max(abs(eta_t));
%             end

%             xi = x;
% %             eta_x  = ifft(1i*kx.*FFTeta);
%             hphi = fft(phiS)/nx.*(k<k_cut); % phi = "phiS" 
%     
% %             FFTeta = fft(eta).*(k<k_cut);
% %             f =  xi + 2i*fft(conj(FFTeta).*(kx>0))/nx; % deep water
% %             f =  xi + 2i*sum( conj(FFTeta.'/nx).*exp(-1i*kx.'.*xi) .*(kx.'>0),2); % deep water
% %             max(abs(imag(f) - eta))
% 
% %             df =  1 + 2i*sum(-1i*kx.'.* conj(FFTeta.'/nx).*exp(-1i*kx.'.*xi) .*(kx.'>0),2); 
% %             df =  1 + 2i*fft(-1i*kx.*conj(FFTeta/nx) .*(kx>0)); 
%             df =  1 + 2*fft(kx.*conj(fft(eta)).*(kx>0&k<k_cut)/nx); 
%             JInv = abs(df).^(-2);
%                         
% %             max(abs( real(fft(2*conj(hphi).*(kx>0))) - phiS )) % test real(omega) = phiS
%             U = conj(-2i*fft(kx.*conj(hphi).*(kx>0)));
%             
%             hb = fft(JInv.*imag(U))/nx;
% %             tf_ = sum( 2i*hb'.*exp(-1i*kx.'.*xi).*(kx.'>0),2) + 1i*hb(1);
%             tf =  fft(2i*conj(hb).*(kx>0&k<k_cut)) + 1i*hb(1);
% %             max(abs(  imag(tf)-JInv.*imag(U)  )) % test imag(tf) = JInv.*imag(U)
%             
%             eta_t = real(df).*imag(tf) + imag(df).*real(tf);
%             phiS_t = real(tf).*real(U) - .5*JInv.*real(U.^2) - g*eta; % phi_t = "phiS_t" 