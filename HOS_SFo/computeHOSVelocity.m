function [u,v,tUsed]=computeHOSVelocity(phiS,eta,x,t,h,xOut,zOut,tOut,M,P)

%Set dimensions of xOut (column) and zOut (row)
if size(xOut,1)==1
    xOut=xOut.';
end
if size(zOut,2)==1
    zOut=zOut.';
end

%Compute wave number vector
Nx=length(x);
Lx=length(x)*(x(2)-x(1));
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
v=zeros(length(xOut),length(zOut),length(tUsed));

%Main time loop
zMat=repmat(zOut,length(xOut),1);
kTab=repmat(k.',1,length(zOut),length(xOut));
xTab=repmat(shiftdim(xOut,-2),length(k),length(zOut),1);
for iT=1:length(tUsed)
    fprintf('t=%8.3f\n',tUsed(iT));
    Axi0=getXi0Modes(phiS(indT(iT),:),eta(indT(iT),:),k,M,h,P);
    indKeep=log10(abs(Axi0)/max(abs(Axi0)))>-8;
%     kTab=repmat(k(indKeep).',1,length(zOut),length(xOut));
%     xTab=repmat(shiftdim(xOut,-2),length(k(indKeep)),length(zOut),1);
    uAmp=repmat(Axi0(indKeep).',1,length(zOut),length(xOut)).*repmat(Ch(indKeep,:),1,1,length(xOut)).*(1i*kTab(indKeep,:,:).*exp(1i*xTab(indKeep,:,:).*kTab(indKeep,:,:)));
    vAmp=repmat(Axi0(indKeep).',1,length(zOut),length(xOut)).*repmat(Sh(indKeep,:),1,1,length(xOut)).*(kTab(indKeep,:,:).*exp(1i*xTab(indKeep,:,:).*kTab(indKeep,:,:)));
    u2=squeeze(real(sum(uAmp))).'/(Nx/2);
    v2=squeeze(real(sum(vAmp))).'/(Nx/2);
    etaMask=(zMat>repmat(interp1(x,eta(indT(iT),:),xOut),1,length(zOut)));
    u2(etaMask)=NaN;
    v2(etaMask)=NaN;
    u(:,:,iT)=u2;
    v(:,:,iT)=v2;
end

end