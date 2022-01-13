function a=getFunctionFromModeAmplitudes(A)

%A   : mxN  double,   N-Fourier amplitudes of the m functions
%a   : mx2N double,   matrix containing the 2N values of m functions

[m,N2]=size(A);
a=ifft([A,zeros(m,1),fliplr(conj(A(:,2:end)))],2*N2,2);

end


