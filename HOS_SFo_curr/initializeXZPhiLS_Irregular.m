function [deltax,deltaz,phiLS]=initializeXZPhiLS_Irregular(Nx,h,k,psi,seed)


%Use Matlab's Mersenne Twister algorithm to generate random phases
rng(double(seed),'twister');
phi=2*pi*rand(size(psi));
complexRand=exp(-1i*phi);

%Compute deltax and deltaz 
dk=k(2)-k(1);
amp=sqrt(2*psi*dk);
fourCoef=Nx*[0;flipud(conj(amp(2:end).*complexRand(2:end)));amp.*complexRand]/2;
deltax=real(ifft(ifftshift(1i*fourCoef)));
deltaz=real(ifft(ifftshift(fourCoef)));

%Compute phiLS
g=9.81;
om=sqrt(k*g.*tanh(k*h));
HphiLS=-1i*om./k;
HphiLS(k==0)=0;
fourCoefPhiLS=Nx*[0;flipud(-HphiLS(2:end).*conj(amp(2:end).*complexRand(2:end)));HphiLS.*amp.*complexRand]/2;
phiLS=real(ifft(ifftshift(fourCoefPhiLS)));


% dx=Lx/Nx;
% x=dx*(0:(Nx-1));
% figure;
% subplot(2,1,1)
% plot(x,eta);grid on
% subplot(2,1,2)
% plot(x,phiS);grid on

end