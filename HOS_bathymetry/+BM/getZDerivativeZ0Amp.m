function dAdz=getZDerivativeZ0Amp(A,n,k,TH)

if mod(n,2)==0
    H=k.^n;
else
    H=(k.^n).*TH;
end
dAdz=A.*H;

end
