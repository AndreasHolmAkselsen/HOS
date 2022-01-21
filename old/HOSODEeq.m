function [phiS_t,eta_t] = HOSODEeq(t,phiS,eta)
    global Tramp nRamp
    wNl = max(0,1-exp(-(t/Tramp)^nRamp));
    g = 9.81;
    [W_lin,W_nl,phiS_x,eta_x] = phiComponentsHOS(phiS,eta);
    W = W_lin+W_nl;
    eta_t  =   W_lin + wNl.*(  W_nl          + eta_x.^2.*W - phiS_x.*eta_x );
    phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*W.^2 );
    
%     eta_t  = -wNl*phiS_x.*eta_x + (1+wNl*eta_x.^2).*W;
%     phiS_t = -.5*wNl*phiS_x.^2 - g*eta + wNl*.5*(1+eta_x.^2).*W.^2;
end
