function [W_lin,W_nl] = phiComponentsHOS(phiS,eta)
    global DO_PADDING taylor
    assert(iscolumn(phiS) && iscolumn(eta));
    M = taylor.M;
    
    % indices:
    % j: spatial/modal index, n: Stokes expansion index, i: z-derivative order (array locatin i+1).
    
    % Allocate memory
    N = size(phiS,1);
    if DO_PADDING
        Nd = N*(M+2)/2; % SFo uses nx*(M+2)/2 instead of nx*(M+1)/2
        phiS = real(ifft(fftPad(fft(phiS),Nd)));
        eta = real(ifft(fftPad(fft(eta),Nd)));
    else
        Nd = N;
    end
    W_jn = zeros(Nd,M);
    hphi_jn = zeros(Nd,M);
    phi_jni = zeros(Nd,M,M+1); % [j,n,i+1] where i is the i'th derivative in z.
    phi_jni(:,1) = phiS;
    
    k = [0:ceil(Nd/2)-1, floor(Nd/2):-1:1]';
    H_ji = k.^(0:M); % [j,i+1] where i is the i'th derivative in z.
    if isfinite(taylor.h), H_ji(:,2:2:M+1) = H_ji(:,2:2:M+1).*tanh(k*taylor.h); end

    for n = 1:M
        % Compute phi^(n)
        for i = 1:n-1
            phi_jni(:,n,1) = phi_jni(:,n,1) - eta.^i.*phi_jni(:,n-i,i+1)/factorial(i);
        end
        % compute new derivatives
        hphi_jn(:,n) = fft(phi_jni(:,n,1));% NB
        for i = 1:(M-n+1)
            phi_jni(:,n,i+1) = real(ifft(fftPad( H_ji(:,i+1).*hphi_jn(:,n),Nd))); % ifft(k.^i.*hphi_jn(:,n));
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
        W_lin = ifft(fftPad(fft(W_lin),N));
        W_nl = ifft(fftPad(fft(W_nl),N));
    end
%     hphi = sum(hphi_jn,2);
%     hphi = [hphi(1:N/2,:);2*hphi(Nd-N/2+1,:);hphi(Nd-N/2+2:Nd,:)];
end