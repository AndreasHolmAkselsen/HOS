function AphiE=computeAPhiE(phiLS,deltaX,deltaZ,k,M,h)

AphiS=getModeAmplitudes(phiLS);
N=length(AphiS);
Nd=(M+2)/2*N;

%Initialize matrices 
phiVol_AA=zeros(M,2*Nd);
phiVolAmp_AA=zeros(M,Nd);
phiVolAmp_AA(1,1:N)=AphiS/N*Nd;
k_AA=(k(2)-k(1))*(0:(Nd-1));
T_AA=tanh(k_AA*h);

%Compute powers of eta_AA with full anti-aliasing treatment
deltaXPow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(deltaX,Nd),M-1,1)]);
deltaZPow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(deltaZ,Nd),M-1,1)]);

%Main loop over HOS-order m
for m=1:M
    for p=0:(m-1)
        for q=0:(m-p)
            nx=m-p-q;
            nz=q;
            phiXZDer_AA=getXZDerivative(phiVolAmp_AA(p+1,:),nx,nz,k_AA,T_AA);
            if m<M
                phiVol_AA(m+1,:)=phiVol_AA(m+1,:)-phiXZDer_AA.*deltaXPow_AA(nx+1,:).*deltaZPow_AA(nz+1,:)/gamma(nx+1)/gamma(nz+1);
            end
        end
    end
    if m<M
        phiVolAmp_AA(m+1,:)=getModeAmplitudes(phiVol_AA(m+1,:));
    end
end

%Return the linear and non linear parts of U and W with N modes only 
AphiE=N/Nd*sum(phiVolAmp_AA(:,1:N));

end

