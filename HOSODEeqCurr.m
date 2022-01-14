function [phiS_t,eta_t] = HOSODEeqCurr(t,phiS,eta)
    global Tramp nRamp x df
    wNl = max(0,1-exp(-(t/Tramp)^nRamp));
    wCurr = wNl;
    g = 9.81;
    [W_lin,W_nl,phiS_x,eta_x] = phiComponentsHOS(phiS,eta);
    W = W_lin+W_nl;
    dfS = df(x + 1i*eta);
    Phi_x =  real(dfS);
    Phi_z = -imag(dfS);
    eta_t  =   W_lin + wNl.*(  W_nl          + eta_x.^2.*W - phiS_x.*eta_x ) - wCurr.*(Phi_x.*eta_x-Phi_z);
    phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*W.^2 )...
             - wCurr.*( Phi_x.*phiS_x + .5*Phi_x.^2 + .5*Phi_z.^2 );
%     phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*W.^2 )...
%              - wCurr.*( Phi_x.*(phiS_x+.5*Phi_x-wNl*eta_x.*W) + Phi_z.*(.5*Phi_z + W) );
    
end
