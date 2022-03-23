function [phiS_t,eta_t] = HOSODEeqCurr(t,phiS,eta)
    global nonLinRamp x df k_cut
    wNl = nonLinRamp(t); 
%     wCurr = wNl;
    wCurr = 1;%wNl;
    g = 9.81;
    
%     % temp!
%     nx = length(eta);
%     dk = 2*pi/(nx*(x(2)-x(1)));
%     if mod(nx,2)==0
%         kx = [0:nx/2-1,-nx/2:-1]'*dk;
%     else
%         kx = [0:(nx-1)/2, -(nx-1)/2:-1]'*dk;
%     end
%     k = abs(kx);
%     k_cut_ = k_cut;
% %     k_cut = inf;
    
    
    [W_lin,W_nl,phiS_x,eta_x] = phiComponentsHOS(phiS,eta);
    W = W_lin+W_nl;
    dfS = df(x + 1i*eta);
    Phi_x =  real(dfS);
    Phi_z = -imag(dfS);
    eta_t  =   W_lin + wNl.*(  W_nl          + eta_x.^2.*W - phiS_x.*eta_x ) - wCurr.*(Phi_x.*eta_x-Phi_z);
    phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*W.^2 )...
             - wCurr.*( Phi_x.*phiS_x + .5*Phi_x.^2 + .5*Phi_z.^2 );    
         
%     heta_t = fft(eta_t); hphiS_t = fft(phiS_t);     
%     eta_t  = ifft(heta_t.*(k>0 & k<k_cut_ & abs(heta_t)>1e-6));
%     phiS_t = ifft(hphiS_t.*(k>0 & k<k_cut_ & abs(hphiS_t)>1e-6));     
         
end
