function [W_lin,W_nl,phiS_x,eta_x,hphi,kx] = phiComponentsHOS(phiS,eta,h,M)
    global x k_cut 
    assert(iscolumn(phiS) && iscolumn(eta));
    
    % indices:
    % j: spatial/modal index, n: Stokes expansion index, i: z-derivative order (array locatin i+1).
    
    kx = getKx(x);
    k = abs(kx);
    
    % Allocate memory
    nx = size(phiS,1);
    W_jn = zeros(nx,M);
    phi_jni = zeros(nx,M,M+1); % [j,n,i+1] where i is the i'th derivative in z.
    phi_jni(:,1) = phiS;
    hphi_jn = zeros(nx,M);

    H_ji = k.^(0:M); % [j,i+1] where i is the i'th derivative in z.
    if isfinite(h), H_ji(:,2:2:M+1) = H_ji(:,2:2:M+1).*tanh(k*h); end
    
    for n = 1:M
        % Compute phi^(n)
        for i = 1:n-1
            phi_jni(:,n,1) = phi_jni(:,n,1) - eta.^i.*phi_jni(:,n-i,i+1)/factorial(i);
        end
        % compute new derivatives
        hphi_jn(:,n) = fft(phi_jni(:,n,1)).*(k<k_cut);% NB
        for i = 1:(M-n+1)
            phi_jni(:,n,i+1) = ifft( H_ji(:,i+1).*hphi_jn(:,n));% ifft(k.^i.*hphi_jn(:,n)); %
        end        
        %compute W^(n)
        for i = 0:n-1
            W_jn(:,n) = W_jn(:,n) + eta.^i.*phi_jni(:,n-i,i+2)/factorial(i);
        end       
    end
    W_jn   = setReal(W_jn,'W');   
    W_lin  = W_jn(:,1);
    W_nl   = sum(W_jn(:,2:M),2);
    eta_x  = setReal( ifft(1i*kx.*fft(eta).*(k<k_cut)), 'eta_x' );
    phiS_x = setReal( ifft(1i*kx.*hphi_jn(:,1)), 'phiS_x' );
    hphi   = sum(hphi_jn,2);        
end



