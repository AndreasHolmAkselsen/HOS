function a=getXZDerivative(A,nx,nz,k,T)

if mod(nz,2)==0
    H=(1i^nx)*k.^(nx+nz);
else
    H=(1i^nx)*(k.^(nx+nz)).*T;
end
dAdz=A.*H;
a=ifft([dAdz,0,fliplr(conj(dAdz(2:end)))]);

end
