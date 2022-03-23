function [U,eta_x,homega,kxp] = phiComponentsConformal(phiS,eta)
    global x k_cut
    assert(iscolumn(phiS) && iscolumn(eta));
    
    nx = size(phiS,1);
    dk = 2*pi/(nx*(x(2)-x(1)));
    if mod(nx,2)==0
        kx = [0:nx/2-1,-nx/2:-1]'*dk;
    else
        kx = [0:(nx-1)/2, -(nx-1)/2:-1]'*dk;
    end
    k = abs(kx);  
    hphiS = fft(phiS)/nx.*(k<k_cut);
    
    % efficient:
    dw = -2i*fft(kx.*conj(hphiS).*(kx>0));
    eta_x  = ifft(1i*kx.*fft(eta).*(k<k_cut));
    U = conj( dw./(1+1i*eta_x) );  
    
    if nargout > 1
        homega = [hphiS(1);  2*conj(hphiS(kx>0))];
        kxp = kx(kx>=0);
    end
    
%     kxp = kx(kx>=0 & k<k_cut); 
%     homega = [hphiS(1);  2*conj(hphiS(kx>0))];
%     eta_x  = ifft(1i*kx.*fft(eta).*(k<k_cut));
%     U = conj( 1./(1i-eta_x).*sum(kxp.*homega.*exp(-1i*kxp.*x.'),1).' );  

%     homega = [hphiS(1);  2*conj(hphiS(kx>0))];
%     heta = fft(eta).*(k<k_cut)/nx;
%     df_ = 1+1i*sum(1i*kx.*heta.*   exp(1i*kx.*x'),1).';
%     dw_ =      sum(-1i*kxp.*homega.*exp(-1i*kxp.*x'),1).';
%     U = conj( dw./df );  
    
end