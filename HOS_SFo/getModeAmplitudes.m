function A=getModeAmplitudes(a)

%a   : mx2N double,   matrix containing the 2N values of m functions
%A   : mxN  double,   N-Fourier amplitudes of the m functions
N2=size(a,2);
A=fft(a,N2,2);
A=A(:,1:N2/2);

end

