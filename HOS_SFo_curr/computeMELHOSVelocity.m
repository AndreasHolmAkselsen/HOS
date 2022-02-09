function [u,w,tUsed]=computeMELHOSVelocity(phiLS,deltaX,deltaZ,eta,alpha,t,h,xOut,zOut,tOut,M)

%Set dimensions of xOut (column) and zOut (row)
if size(xOut,1)==1
    xOut=xOut.';
end
if size(zOut,2)==1
    zOut=zOut.';
end

%Compute wave number vector
Nx=length(alpha);
Lx=Nx*(alpha(2)-alpha(1));
dk=2*pi/Lx;
k=dk*(0:(Nx/2-1));

%Precompute cosh and sinh matrices
C=k.'*(zOut+h);
D=k.'*repmat(h,1,length(zOut));
limMask=(C<10)|(D<10);
Ch=exp(C-D);   %Simplified expression to avoid infty/infty terms
Ch(limMask)=cosh(C(limMask))./cosh(D(limMask));
Sh=Ch;
Sh(limMask)=sinh(C(limMask))./cosh(D(limMask));

%Find time indices of interest
indT=find(ismembertol(t,tOut,1e-10));
tUsed=t(indT);

%Define size of velocity outputs
u=zeros(length(xOut),length(zOut),length(tUsed));
w=zeros(length(xOut),length(zOut),length(tUsed));

%Main time loop
zMat=repmat(zOut,length(xOut),1);
kTab=repmat(k.',1,length(zOut),length(xOut));
xTab=repmat(shiftdim(xOut,-2),length(k),length(zOut),1);
for iT=1:length(tUsed)
%     fprintf('t=%8.3f\n',tUsed(iT));
    deltaXPrime=mean(deltaX(indT(iT),:));
    AphiE=computeAPhiE(phiLS(indT(iT),:),deltaX(indT(iT),:)-deltaXPrime,deltaZ(indT(iT),:),k,M,h);
    indKKeep=(log10(abs(AphiE)/max(abs(AphiE)))>-8);
    uAmp=repmat(AphiE(indKKeep).',1,length(zOut),length(xOut)).*repmat(Ch(indKKeep,:),1,1,length(xOut)).*...
        (1i*kTab(indKKeep,:,:).*exp(1i*(xTab(indKKeep,:,:)-deltaXPrime).*kTab(indKKeep,:,:)));
    vAmp=repmat(AphiE(indKKeep).',1,length(zOut),length(xOut)).*repmat(Sh(indKKeep,:),1,1,length(xOut)).*...
        (kTab(indKKeep,:,:).*exp(1i*(xTab(indKKeep,:,:)-deltaXPrime).*kTab(indKKeep,:,:)));
    u2=squeeze(real(sum(uAmp))).'/(Nx/2);
    w2=squeeze(real(sum(vAmp))).'/(Nx/2);
    etaMask=(zMat>repmat(interp1(alpha,eta(indT(iT),:),xOut),1,length(zOut)));
    u2(etaMask)=NaN;
    w2(etaMask)=NaN;
    u(:,:,iT)=u2;
    w(:,:,iT)=w2;
end

end

