function [phiS_t,eta_t] = HOSODEeqCurr(t,phiS,eta)
    global nonLinRamp x df
    wNl = nonLinRamp(t); 
    wCurr = 1;%wNl;
    g = 9.81;
    [W_lin,W_nl,phiS_x,eta_x] = phiComponentsHOS(phiS,eta);
    W = W_lin+W_nl;
    dfS = df(x + 1i*eta);
    Phi_x =  real(dfS);
    Phi_z = -imag(dfS);
    eta_t  =   W_lin + wNl.*(  W_nl          + eta_x.^2.*W - phiS_x.*eta_x ) - wCurr.*(Phi_x.*eta_x-Phi_z);
    phiS_t = - g*eta + wNl.*( -.5*phiS_x.^2  + .5*(1+eta_x.^2).*W.^2 )...
             - wCurr.*( Phi_x.*phiS_x + .5*Phi_x.^2 + .5*Phi_z.^2 );    
end
