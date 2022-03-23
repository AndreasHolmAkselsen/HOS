function [U,eta_x] = phiComponentsConformalFiniteDepth(phiS,eta,eta_t)
    global x k_cut 
    assert(iscolumn(phiS) && iscolumn(eta));

    kx = getKx(x);
    k = abs(kx);  
    hphiS = fft(phiS)/nx.*(k<k_cut);
    
    
    
%     hetap = heta(kx>0).';
%     kp = kx(kx>0).';
%     xi = x; %nb!
% %     f =  zeta + 2i*sum( conj(hetap).*exp(1i*kp.*zeta) ,2); % deep water
% %     f =@(zeta) zeta - 1i*sum( heta.*exp(1i*k.*(zeta+1i*H))./sinh(k*H) ,3);
%     df =  1 + 2i*sum(1i*kp.*conj(hetap).*exp(1i*kp.*ix) ,2);
%     
    % efficient:
    U = conj(-2i*fft(kx.*conj(hphiS).*(kx>0)));
    eta_x  = ifft(1i*kx.*fft(eta).*(k<k_cut));   

    
end