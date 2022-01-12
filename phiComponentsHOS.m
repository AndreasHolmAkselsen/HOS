function [W,phiS_x,eta_x,hphi,kx] = phiComponentsHOS(phiS,eta)
    global N dx k_cut
    
    assert(iscolumn(phiS) && iscolumn(eta));

    nx = size(phiS,1);
    
    W_jn = zeros(nx,N);
    
    % phi_jni in [j,n,i+1] where i is the i'th derivative in z.
    phi_jni = zeros(nx,N,N+1);
    phi_jni(:,1) = phiS;
    
    hphi_jn = zeros(nx,N);
    
    dk = 2*pi/(nx*dx);
    if mod(nx,2)==0
        kx = [0:nx/2-1,-nx/2:-1]'*dk;
    else
        kx = [0:(nx-1)/2, -(nx-1)/2:-1]'*dk;
    end
    k = abs(kx);
    
    for n = 1:N
        % Compute phi^(n)
        for i = 1:n-1
            assert(any(phi_jni(:,n-i,i)~=0))
            phi_jni(:,n,1) = phi_jni(:,n,1) - eta.^i.*phi_jni(:,n-i,i+1)/factorial(i);
        end
        
        % compute new derivatives
        hphi_jn(:,n) = fft(phi_jni(:,n,1)).*(k<k_cut);% NB
        for i = 1:(N-n+1)
            phi_jni(:,n,i+1) = ifft(k.^i.*hphi_jn(:,n));
        end        
        
        %compute W^(n)
        for i = 0:n-1
            assert(any(phi_jni(:,n-i,i+1)~=0))
            W_jn(:,n) = W_jn(:,n) + eta.^i.*phi_jni(:,n-i,i+2)/factorial(i);
        end       
    end
        
    eta_x  = setReal( ifft(1i*kx.*fft(eta).*(k<k_cut)), 'eta_x' );
    phiS_x = setReal( ifft(1i*kx.*hphi_jn(:,1)), 'phiS_x' );
    W      = setReal( sum(W_jn,2), 'W' );   
    hphi   = sum(hphi_jn,2);
end