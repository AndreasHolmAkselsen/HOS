function dAdz=getZAddDerivativeAmp(A,n,k,Nder)

if mod(n,4)==1||mod(n,4)==2
    H=-k.^n;
else
    H=k.^n;
end
dAdz=A.*H;
dAdz(Nder:end)=0;

end
