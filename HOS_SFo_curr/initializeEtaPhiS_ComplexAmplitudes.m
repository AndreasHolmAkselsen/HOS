function [eta,phiS]=initializeEtaPhiS_ComplexAmplitudes(k,compAmp,h)

%Compute deltax and deltaz 
fourCoefz=length(k)*[0;flipud(conj(compAmp(2:end)));compAmp];
eta=real(ifft(ifftshift(fourCoefz)));

%Compute phiLS
g=9.81;
om=sqrt(k*g.*tanh(k*h));
HphiS=-1i*g./om;
% HphiS=-1i*om./k;
HphiS(k==0)=0;
fourCoefPhiS=length(k)*[0;flipud(-HphiS(2:end).*conj(compAmp(2:end)));HphiS.*compAmp];
phiS=real(ifft(ifftshift(fourCoefPhiS)));

end