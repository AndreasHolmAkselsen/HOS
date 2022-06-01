function A=getCosModeAmplitudes(a)

%Decompose a function f(x) on [0,L] into cos(pi*k*x/L) modes
%   with amplitudes a(k) for k=0..N, where N is assumed to be even
%
%a   : m x (N+1) double,   matrix containing the N+1 values of m functions
%A   : m x (N+1) double,   amplitudes of the N+1 cos-modes of the m functions

N=size(a,2)-1;
aPer=[a,fliplr(a(:,2:end-1))];
A=fft(aPer,size(aPer,2),2);
A=real(A(:,1:N+1));

end

