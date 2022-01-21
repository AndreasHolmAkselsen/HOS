function [eta,hphi,kx] = HOSODEeqCurr_linearSteadyState(x,df,nIt)
g = 9.81;

dfS = df(x);
Phi_x =  real(dfS);
Phi_z = -imag(dfS);
assert(all(Phi_z==0))
nx = length(x);
L = nx*(x(2)-x(1));
dk = 2*pi/L;
if mod(nx,2)==0
    kx = [0:nx/2-1,-nx/2:-1]'*dk;
else
    kx = [0:(nx-1)/2, -(nx-1)/2:-1]'*dk;
end
k = abs(kx);

phi_x = 0;
etaOld = 0;
for it = 1:nIt
    eta = -.5/g*Phi_x.^2 - 1/g*phi_x.*Phi_x;
    eta_x = ifft(1i*kx.*fft(eta));
    hphi = fft(eta_x.*Phi_x)./k ; hphi(k==0) = 0;
    phi_x = ifft(1i*kx.*hphi);
    
    fprintf('iteration %d: maxChange = %g\n',it, max(abs(etaOld-eta)));
    if max(abs(etaOld-eta)) < 5e-5,  break; end
    etaOld = eta;
end

    


end