function [W_lin,W_nl] = phiComponentsHOS(phiS,eta)
    global M H DO_PADDING kx
    assert(iscolumn(phiS) && iscolumn(eta));
    
    % indices:
    % j: spatial/modal index, n: Stokes expansion index, i: z-derivative order (array locatin i+1).
    
    % Allocate memory
    N = size(phiS,1);
    if DO_PADDING
        Nd = N*(M+2)/2; % SFo uses nx*(M+2)/2 instead of nx*(M+1)/2
        phiS = real(ifftPad(fft(phiS),Nd));
        eta = real(ifftPad(fft(eta),Nd));
    else
        Nd = N;
    end
    W_jn = zeros(Nd,M);
    hphi_jn = zeros(Nd,M);
    phi_jni = zeros(Nd,M,M+1); % [j,n,i+1] where i is the i'th derivative in z.
    phi_jni(:,1) = phiS;
    
    
%     kx = getKx(x);
%     dk = 2*pi/(N*(x(2)-x(1)));
    if mod(Nd,2)==0
        kPad = [0:Nd/2-1,Nd/2:-1:1]'*kx(2);
    else
        kPad = [0:(Nd-1)/2, (Nd-1)/2:-1:1]'*kx(2);
    end

    H_ji = kPad.^(0:M); % [j,i+1] where i is the i'th derivative in z.
    if isfinite(H), H_ji(:,2:2:M+1) = H_ji(:,2:2:M+1).*tanh(kPad*H); end

    for n = 1:M
        % Compute phi^(n)
        for i = 1:n-1
            phi_jni(:,n,1) = phi_jni(:,n,1) - eta.^i.*phi_jni(:,n-i,i+1)/factorial(i);
        end
        % compute new derivatives
        hphi_jn(:,n) = fft(phi_jni(:,n,1));% NB
        for i = 1:(M-n+1)
            phi_jni(:,n,i+1) = real(ifftPad( H_ji(:,i+1).*hphi_jn(:,n),Nd)); % ifft(k.^i.*hphi_jn(:,n));
        end        
        %compute W^(n)
        for i = 0:n-1
            W_jn(:,n) = W_jn(:,n) + eta.^i.*phi_jni(:,n-i,i+2)/factorial(i);
        end       
    end
%     W_jn   = real(W_jn);   
    W_lin  = W_jn(:,1);
    W_nl   = sum(W_jn(:,2:M),2);
    if DO_PADDING
        W_lin = ifftPad(fft(W_lin),N);
        W_nl = ifftPad(fft(W_nl),N);
    end
%     hphi = sum(hphi_jn,2);
%     hphi = [hphi(1:N/2,:);2*hphi(Nd-N/2+1,:);hphi(Nd-N/2+2:Nd,:)];
end