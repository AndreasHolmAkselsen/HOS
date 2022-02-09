function psi=computeWaveSpectralDensity_k_JONSWAP(k,Hs,Tp,gam,h)

g=9.81;
f=sqrt(k*g.*tanh(k*h))/(2*pi);
S=level1.computeWaveSpectralDensity_JONSWAP(f,Hs,Tp,gam);
jaco=g/(8*pi^2)./f.*(tanh(k*h)+k*h.*sech(k*h));
psi=S.*jaco;
psi(k==0)=0;
psi=psi*(Hs/(4*sqrt(sum(psi)*(k(2)-k(1)))))^2;

end

