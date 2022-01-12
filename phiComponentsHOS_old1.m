function [W,phiS_x,eta_x,phi,k] = phiComponentsHOS(phiS,eta,N,dx)
    
    assert(iscolumn(phiS) && iscolumn(eta));

    
    nx = size(phiS,1);
    
    [W_jn,phi_jn] = deal(zeros(nx,N));
    phi_jn(:,1) = phiS;
    hphi_jn = zeros(nx,N-1);
    phi_z_jni = nan(nx,N-1,N);
    hphi_jn(:,1) = fft(phiS);
    
    
    dk = 1/(nx*dx);
	k = [0:nx/2-1,-nx/2:-1]'*dk;
    for n = 2:N
        % compute new derivatives
        for i = 1:(N-n+1)
            phi_z_jni(:,n-1,i) = ifft(abs(k).^i.*hphi_jn(:,n-1));
        end
        for i = 1:n-1
            assert(~any(isnan(phi_z_jni(:,n-i,i))))
            if nargout > 1, phi_jn(:,n) = phi_jn(:,n) - eta.^i.*phi_z_jni(:,n-i,i)/factorial(i); end
            W_jn(:,n) = W_jn(:,n) + eta.^i.*phi_z_jni(:,n-i,i+1)/factorial(i);
        end
        if n~=N, hphi_jn(:,n) = fft(phi_jn(:,n)); end
    end
    

        
    testReal = @(x) assert( all( abs(imag(x)./(abs(x)+1e-9)) < 1e-6,'all') );
    
    eta_x = ifft(1i*k.*fft(eta));
    phiS_x = ifft(1i*k.*hphi_jn(:,1));
    W = sum(W_jn,2);
    
    % Remove imaginary zeros
	
    testReal(phiS_x); phiS_x = real(phiS_x);
    testReal(eta_x); eta_x = real(eta_x);
    
    if nargout > 1, phi = sum(phi_jn,2); testReal(phi); phi = real(phi); end
    
    
%     hphiS2 = hphi_jn.*exp(abs(k).*eta);
%     figure, plot(k,real(hphiS),k,real(hphiS2),'--',k,imag(hphiS),k,imag(hphiS2),'--')
end