function [phiS_t,eta_t] = HOSODEeq(phiS,eta)

    wNl = 1;
    g = 9.81;
    [W,phiS_x,eta_x] = phiComponentsHOS(phiS,eta);
    eta_t  = -wNl*phiS_x.*eta_x + (1+wNl*eta_x.^2).*W;
    phiS_t = -.5*wNl*phiS_x.^2 - g*eta + wNl*.5*(1+eta_x.^2).*W.^2;
end

