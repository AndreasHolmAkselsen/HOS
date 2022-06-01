function dAdx=getXAddDerivativeX0Amp(A,n,k,T)

if mod(n,2)==0
    H=(-k).^n;
else
    H=((-k).^n).*T.^(sign(n));
end
H(1)=0;
dAdx=A.*H;

end
