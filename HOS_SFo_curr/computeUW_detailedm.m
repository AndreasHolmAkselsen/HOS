function [Ulin,Unl,Wlin,Wnl,AphiE]=computeUW_Detailed(deltaX,deltaZ,AphiS,k,M,h)

N=length(AphiS);
Nd=(M+2)/2*N;

%Initialize matrices 
phiVol_AA=zeros(M,2*Nd);
phiVolAmp_AA=zeros(M,Nd);
U_AA=zeros(M,2*Nd);
W_AA=zeros(M,2*Nd);
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
            if nx>=1
                U_AA(m,:)=U_AA(m,:)+phiXZDer_AA.*deltaXPow_AA(nx,:).*deltaZPow_AA(nz+1,:)/gamma(nx)/gamma(nz+1);
            end
            if nz>=1
                W_AA(m,:)=W_AA(m,:)+phiXZDer_AA.*deltaXPow_AA(nx+1,:).*deltaZPow_AA(nz,:)/gamma(nx+1)/gamma(nz);
            end
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
AU_AA=getModeAmplitudes(U_AA);
AW_AA=getModeAmplitudes(W_AA);
U=N/Nd*getFunctionFromModeAmplitudes(AU_AA(:,1:N));
W=N/Nd*getFunctionFromModeAmplitudes(AW_AA(:,1:N));
phiVol=N/Nd*getFunctionFromModeAmplitudes(phiVolAmp_AA(:,1:N));
Ulin=U(1,:);
Unl=sum(U(2:end,:),1);
Wlin=W(1,:);
Wnl=sum(W(2:end,:),1);
AphiE=getModeAmplitudes(sum(phiVol,1));

end
