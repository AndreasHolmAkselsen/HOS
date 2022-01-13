function [axConv,azConv]=computeAccConv(deltaX,deltaZ,AphiE,k,M,h,abcd)

N=length(AphiE);
Nd=(M+2)/2*N;

%Initialize matrices 
ax_AA=zeros(M-1,2*Nd);
az_AA=zeros(M-1,2*Nd);
phiEAmp_AA=zeros(1,Nd);
phiEAmp_AA(1,1:N)=AphiE/N*Nd;
k_AA=(k(2)-k(1))*(0:(Nd-1));
T_AA=tanh(k_AA*h);

%Precompute powers of deltaX and deltaZ with full anti-aliasing treatment
deltaXPow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(deltaX,Nd),M-1,1)]);
deltaZPow_AA=cumprod([ones(1,2*Nd);repmat(zeroPadding(deltaZ,Nd),M-1,1)]);

%Precompute x- and z-derivatives of the Eulerian potential at z=0
phiDotXZDer_AA=zeros(M*(M+1)-1,2*Nd);
for b=0:M
    for a=0:M
        if (a~=0)||(b~=0)
            phiDotXZDer_AA(b*(M+1)+a,:)=getXZDerivative(phiEAmp_AA,a,b,k_AA,T_AA);
%             fprintf('%3d --> %3d %3d\n',b*(M+1)+a,a,b);
        end
    end
end
% fprintf('\n\n');

%Compute Taylor expansion
for n=0:M-2
    for p=0:n
        dphiProdX=zeros(1,2*Nd);
        dphiProdZ=zeros(1,2*Nd);
        for i0=1:size(abcd{n+1,p+1},1)
            a=abcd{n+1,p+1}(i0,1);
            b=abcd{n+1,p+1}(i0,2);
            c=abcd{n+1,p+1}(i0,3);
            d=abcd{n+1,p+1}(i0,4);
            w=abcd{n+1,p+1}(i0,5);
            dphiProdX=dphiProdX+w*(phiDotXZDer_AA(b*(M+1)+(a+1),:).*phiDotXZDer_AA(d*(M+1)+(c+2),:)+...
                phiDotXZDer_AA((b+1)*(M+1)+a,:).*phiDotXZDer_AA((d+1)*(M+1)+(c+1),:));
            dphiProdZ=dphiProdZ+w*(phiDotXZDer_AA(b*(M+1)+(a+1),:).*phiDotXZDer_AA((d+1)*(M+1)+(c+1),:)+...
                phiDotXZDer_AA((b+1)*(M+1)+a,:).*phiDotXZDer_AA((d+2)*(M+1)+c,:));
            %             fprintf('%3d %3d --> %3d x (%d,%d) (%d,%d) + (%d,%d) (%d,%d)\n',p,n-p,w,a+1,b,c+2,d,a,b+1,c+1,d+1),
        end
        %         fprintf('\n');
        ax_AA(n+1,:)=ax_AA(n+1,:)+dphiProdX.*deltaXPow_AA(p+1,:).*deltaZPow_AA(n-p+1,:)/gamma(p+1)/gamma(n-p+1);
        az_AA(n+1,:)=az_AA(n+1,:)+dphiProdZ.*deltaXPow_AA(p+1,:).*deltaZPow_AA(n-p+1,:)/gamma(p+1)/gamma(n-p+1);
    end
%     fprintf('\n');
end
% keyboard

%Return the linear and non linear parts of U and W with N modes only 
Aax_AA=getModeAmplitudes(ax_AA);
axConv=N/Nd*getFunctionFromModeAmplitudes(sum(Aax_AA(:,1:N),1));
Aaz_AA=getModeAmplitudes(az_AA);
azConv=N/Nd*getFunctionFromModeAmplitudes(sum(Aaz_AA(:,1:N),1));

end
