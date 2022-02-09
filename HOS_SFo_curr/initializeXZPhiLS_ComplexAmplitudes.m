function [deltax,deltaz,phiLS,dphiLSdt]=initializeXZPhiLS_ComplexAmplitudes(k,compAmp,h)

%Compute deltax and deltaz 
fourCoefz=length(k)*[0;flipud(conj(compAmp(2:end)));compAmp];
deltaz=real(ifft(ifftshift(fourCoefz)));
fourCoefx=length(k)*[0;flipud(conj(1i*compAmp(2:end)));1i*compAmp];
deltax=real(ifft(ifftshift(fourCoefx)));

%Compute phiLS
g=9.81;
om=sqrt(k*g.*tanh(k*h));
HphiLS=-1i*om./k;
HphiLS(k==0)=0;
fourCoefPhiLS=length(k)*[0;flipud(-HphiLS(2:end).*conj(compAmp(2:end)));HphiLS.*compAmp];
phiLS=real(ifft(ifftshift(fourCoefPhiLS)));
dphiLSdt=-g*deltaz;

% dx=Lx/Nx;
% x=dx*(0:(Nx-1));
% figure;
% subplot(2,1,1)
% plot(x,eta);grid on
% subplot(2,1,2)
% plot(x,phiS);grid on

end