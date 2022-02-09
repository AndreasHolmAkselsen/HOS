function xi2=H2Operator(xi,eta,eta2,k,M,h)

Axi=getModeAmplitudes(xi);
N=length(Axi);
Nd=(M+1)/2*N;  %Full anti-aliasing

%Initialize matrices 
xiVol_AA=zeros(M,2*Nd);
xiVolAmp_AA=zeros(M,Nd);
xiVolAmp_AA(1,1:N)=Axi/N*Nd;
xi2_AA=getFunctionFromModeAmplitudes(xiVolAmp_AA(1,:));
k_AA=(k(2)-k(1))*(0:(Nd-1));
T_AA=tanh(k_AA*h);

%Compute powers of eta and eta2 with full anti-aliasing treatment
etaPow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(eta,Nd),M-1,1)]);
eta2Pow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(eta2,Nd),M-1,1)]);

%Main loop over HOS-order m
for m=2:M
    for p=1:(m-1)
        xiZDer_AA=getZDerivative(xiVolAmp_AA(m-p,:),p,k_AA,T_AA);
        xiVol_AA(m,:)=xiVol_AA(m,:)-xiZDer_AA.*etaPow_AA(p+1,:)/gamma(p+1);
        xi2_AA=xi2_AA+xiZDer_AA.*eta2Pow_AA(p+1,:)/gamma(p+1);
    end
    xi2_AA=xi2_AA+xiVol_AA(m,:);
    xiVolAmp_AA(m,:)=getModeAmplitudes(xiVol_AA(m,:));
end

%Return xi2 with N modes only 
Axi2_AA=getModeAmplitudes(xi2_AA);
xi2=N/Nd*getFunctionFromModeAmplitudes(Axi2_AA(1:N));

end

