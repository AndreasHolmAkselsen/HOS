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

        case 'decayingConformal'  
            nx = numel(x);
            kx = getKx(x); % x->xi
            k = abs(kx);
            if ~isfinite(H), H = realmax; end
%             
%             Lsin = -2./(exp(2*kx.*H)-1-2*(k==0)).*(k<k_cut&k~=0);
            Lsin = -2./(exp(2*kx.*H)-1); Lsin(1) = 0;
            Lcos = 2./(exp(2*kx.*H)+1);
            heta = fft(eta).*(k<k_cut)/nx;
            df =  1 - ifft( kx.*heta.*Lsin)*nx;
            if any(real(df) < 0)% downward mapping -> wave breaking
                [eta_t,phiS_t] = deal(nan(size(k)));  return
            end
            JInv = abs(df).^(-2);
            U = conj(ifft(1i.*kx.*fft(phiS).*(k<k_cut).*Lcos));
            hb = fft(JInv.*imag(U)).*(k<k_cut)/nx;
            f0 =  1i*sum(imag(hb.*conj(heta)).*kx./(sinh(kx*H)+(k==0)).^2); % enforces mean(x_t) = 0
            tf =  1i*f0 + 1i*ifft(hb.*Lsin)*nx;

            eta_t = imag(tf.*df);
            phiS_t = real(tf).*real(U) - .5*JInv.*real(U.^2) - g*eta; 
%             phiS_t = real(conj(U).*tf) - .5*JInv.*abs(U).^2 - g*eta;

        case 'simpleMethod' % the much simpler but erronious method that "curls" with depth.
            [w_lin,~,~,eta_x] = phiComponentsHOS(phiS,eta,H,1);
            U = phiComponentsConformal(phiS,eta);
            w = imag(U); u = real(U);
            w_nl = w-w_lin;
            eta_t  = w_lin +  wNl.*(w_nl - u.*eta_x);
            phiS_t = - g*eta + wNl.*(eta_t.*w -.5*(u.^2+w.^2));
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
