function [W,phiVol]=computeW_detailed(eta,AphiS,k,M,h)

N=length(AphiS);
Nd=(M+2)/2*N;

%Initialize matrices 
phiVol_AA=zeros(M,2*Nd);
phiVolAmp_AA=zeros(M,Nd);
W_AA=zeros(M,2*Nd);
phiVolAmp_AA(1,1:N)=AphiS/N*Nd;
k_AA=(k(2)-k(1))*(0:(Nd-1));
T_AA=tanh(k_AA*h);

%Compute powers of eta_AA with full anti-aliasing treatment
etaPow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(eta,Nd),M-1,1)]);

%Main loop over HOS-order m
for m=1:M
    for p=0:(m-1)
        phiZDer_AA=getZDerivative(phiVolAmp_AA(m-p,:),p+1,k_AA,T_AA);
        W_AA(m,:)=W_AA(m,:)+phiZDer_AA.*etaPow_AA(p+1,:)/gamma(p+1);
        if m<M
            phiVol_AA(m+1,:)=phiVol_AA(m+1,:)-phiZDer_AA.*etaPow_AA(p+2,:)/gamma(p+2);
        end
    end
    if m<M
        phiVolAmp_AA(m+1,:)=getModeAmplitudes(phiVol_AA(m+1,:));
    end
end

%Return W with N modes only 
AW_AA=getModeAmplitudes(W_AA);
W=N/Nd*getFunctionFromModeAmplitudes(AW_AA(:,1:N));
phiVol=N/Nd*getFunctionFromModeAmplitudes(phiVolAmp_AA(:,1:N));

end
