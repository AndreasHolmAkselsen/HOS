function [eta,phiS]=initializeEtaPhiS_Irregular(Nx,h,k,psi,seed)


%Use Matlab's Mersenne Twister algorithm to generate random phases
rng(double(seed),'twister');
phi=2*pi*rand(size(psi));
complexRand=exp(-1i*phi);

%Compute deltax and deltaz 
dk=k(2)-k(1);
amp=sqrt(2*psi*dk);
fourCoef=Nx*[0;flipud(conj(amp(2:end).*complexRand(2:end)));amp.*complexRand]/2;
eta=real(ifft(ifftshift(fourCoef)));

%Compute phiLS
g=9.81;
om=sqrt(k*g.*tanh(k*h));
HphiS=-1i*om./k;
HphiS(k==0)=0;
fourCoefPhiS=Nx*[0;flipud(-HphiS(2:end).*conj(amp(2:end).*complexRand(2:end)));HphiS.*amp.*complexRand]/2;
phiS=real(ifft(ifftshift(fourCoefPhiS)));

end