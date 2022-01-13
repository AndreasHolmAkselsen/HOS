function psi=computeWaveSpectralDensity_k_Gaussian(k,Hs,k0,sigmak)

psi=exp(-(k-k0).^2/(2*sigmak^2));
psi(k==0)=0;
psi(psi<1e-6)=0;
psi=psi*(Hs/(4*sqrt(sum(psi)*(k(2)-k(1)))))^2;

end

