function a=getFunctionFromCosModes(A)

%Reconstruct a function f(x) on [0,L] from cos(pi*k*x/L) modes
%   with amplitudes a(k) for k=0..N, where N is assumed to be even
%
%A   : m x (N+1) double,   amplitudes of the N+1 cos-modes of the m functions
%a   : m x (N+1) double,   matrix containing the N+1 values of m functions

N=size(A,2)-1;
aPer=ifft([A,fliplr(A(:,2:end-1))],2*N,2);
a=aPer(:,1:N+1);

end

