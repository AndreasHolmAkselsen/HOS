function eta=computeEta(deltaX,z,k,M)

Az=getModeAmplitudes(z);
N=length(Az);
Nd=(M+2)/2*N;

%Initialize matrices 
eta_AA=zeros(M,2*Nd);
etaAmp_AA=zeros(M,Nd);
etaAmp_AA(1,1:N)=Az/N*Nd;
k_AA=(k(2)-k(1))*(0:(Nd-1));
T_AA=ones(size(k_AA));

%Compute powers of deltaX with full anti-aliasing treatment
deltaXPow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(deltaX,Nd),M-1,1)]);

%Main loop over HOS-order m
for m=2:M
    for p=1:(m-1)
        etaXDer_AA=getXZDerivative(etaAmp_AA(m-p,:),p,0,k_AA,T_AA);
        eta_AA(m,:)=eta_AA(m,:)-etaXDer_AA.*deltaXPow_AA(p+1,:)/gamma(p+1);
    end
    if m<M
        etaAmp_AA(m,:)=getModeAmplitudes(eta_AA(m,:));
    end
end

%Return eta(m) with N modes only
Aeta_AA=getModeAmplitudes(eta_AA);
eta=N/Nd*getFunctionFromModeAmplitudes(Aeta_AA(:,1:N));
eta(1,:)=z;

end

