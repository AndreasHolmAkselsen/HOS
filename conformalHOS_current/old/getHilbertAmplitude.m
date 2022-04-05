function hilbertAmp=getHilbertAmplitude(x)

    %Compute the amplitude a_H(t) of an amplitude modulated signal by using Hilbert transform, i.e.
    %   ts.value = a_H(t)*cos(omega*t+phi)
    %   where a_H(t) is a slowly varying envelop amplitude 
    %         cos(omaga*t+phi) is a sinusoidal carrier
    %   It is assumed that a_H(t) has no frequency content above the carrier frequency omega (Bedrosian theorem).
    %

    
    %Compute Hilbert transform
%     dom=2*pi/(dt*size(x,1));
%     ommax=pi/dt;

    dom=2/size(x,1);
    ommax=1;

    om=(-ommax:dom:(ommax-dom)).';
    Fts=fftshift(fft(x));
    Hts=real(ifft(ifftshift(-1i*repmat(sign(om),1,size(x,2)).*Fts)));
    hilbertAmp=sqrt(x.^2+Hts.^2);
end