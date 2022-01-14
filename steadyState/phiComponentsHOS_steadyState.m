function [phiS_x,eta_x,phiS,hphi,kx] = phiComponentsHOS_steadyState(W,eta)
    global M dx k_cut
    
    assert(iscolumn(W) && iscolumn(eta));

    nx = size(W,1);
    
    
%     W_jn = zeros(nx,M);
    phiS_jn = zeros(nx,M);
    
    % phi_jni in [j,n,i+1] where i is the i'th derivative in z.
    phi_jni = zeros(nx,M,M+1);
    hphi_jn = zeros(nx,M);
    
    dk = 2*pi/(nx*dx);
    if mod(nx,2)==0
        kx = [0:nx/2-1,-nx/2:-1]'*dk;
    else
        kx = [0:(nx-1)/2, -(nx-1)/2:-1]'*dk;
    end
    k = abs(kx);
    
%     phi_jni(:,1,1) = phiS;
    hphi1 = fft(W)./k; hphi1(k>=k_cut | k==0) = 0;
    phi_jni(:,1,1) = ifft(hphi1);
    
    for n = 1:M
        % Compute phi^(n)
        for i = 1:n-1
%             assert(any(phi_jni(:,n-i,i)~=0)) % test that necessary components have been computed
            phi_jni(:,n,1) = phi_jni(:,n,1) - eta.^i.*phi_jni(:,n-i,i+1)/factorial(i);
        end
        
        % compute new derivatives
        hphi_jn(:,n) = fft(phi_jni(:,n,1)).*(k<k_cut);% MB
        for i = 1:(M-n+1)
            phi_jni(:,n,i+1) = ifft(k.^i.*hphi_jn(:,n));
        end        
        %compute W^(n)
        for i = 0:n-1
%             assert(any(phi_jni(:,n-i,i+1)~=0)) % test that necessary components have been computed
%             W_jn(:,n) = W_jn(:,n) + eta.^i.*phi_jni(:,n-i,i+2)/factorial(i);
            phiS_jn(:,n) = phiS_jn(:,n) + eta.^i.*phi_jni(:,n-i,i+1)/factorial(i);
        end       
    end
%     phi_ji = sum(phi_jni,2);
%     for i = 0:M
%         phiS = phiS + eta.^i.*phi_ji(:,:,i)/factorial(i);
%     end    
    phiS = sum(phiS_jn,2);
    eta_x  = setReal(ifft(1i*kx.*fft(eta).*(k<k_cut)),'eta');
    phiS_x = setReal(ifft(1i*kx.*fft(phiS).*(k<k_cut)),'phiS');
    hphi   = sum(hphi_jn,2);
        
end