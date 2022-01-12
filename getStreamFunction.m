function [phi,psi] = getStreamFunction(dx,z,hphi)

assert(iscolumn(hphi));
assert(isrow(z));

nx = length(hphi);
dk = 2*pi/(nx*dx);
if mod(nx,2)==0
    k = [0:nx/2-1,-nx/2:-1]'*dk;
else
    k = [0:(nx-1)/2, -(nx-1)/2:-1]'*dk;
end

phi = setReal(ifft(hphi.*exp(abs(k).*z),[],1),'phi');
psi = setReal(ifft(1i*sign(k).*hphi.*exp(abs(k).*z),[],1),'psi');

end