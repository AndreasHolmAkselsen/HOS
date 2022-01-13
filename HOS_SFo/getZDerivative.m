function a=getZDerivative(A,n,k,T)

if mod(n,2)==0
    H=k.^n;
else
    H=(k.^n).*T;
end
dAdz=A.*H;
a=ifft([dAdz,0,fliplr(conj(dAdz(2:end)))]);

end
