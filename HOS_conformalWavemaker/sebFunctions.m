function sebFunctions(a)

% A = getCosModeAmplitudes(a);
% a_ = getFunctionFromCosModes(A);

A = getSinModeAmplitudes(a);
a_ = getFunctionFromSinModes(A);

figure('color','w');
n = length(a);
plot(0:n-1,a,0:n-1,a_,'--')


end



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



function A=getSinModeAmplitudes(a)

%Decompose a function f(x) on [0,L] into cos(pi*k*x/L) modes
%   with amplitudes a(k) for k=0..N, where N is assumed to be even
%
%a   : m x (N+1) double,   matrix containing the N+1 values of m functions
%A   : m x (N+1) double,   amplitudes of the N+1 cos-modes of the m functions

N=size(a,2)-1;
% aPer=[a,-fliplr(a(2:end-1))];
aPer=1i*[a,-fliplr(a(2:end-1))]; % NB! correction.
A=fft(aPer,size(aPer,2),2);
A=real(A(:,1:N+1));

end


function a=getFunctionFromSinModes(A)

%Reconstruct a function f(x) on [0,L] from cos(pi*k*x/L) modes
%   with amplitudes a(k) for k=0..N, where N is assumed to be even
%
%A   : m x (N+1) double,   amplitudes of the N+1 cos-modes of the m functions
%a   : m x (N+1) double,   matrix containing the N+1 values of m functions

N=size(A,2)-1;
aPer=ifft(-1i*[A,-fliplr(A(2:end-1))],2*N,2);
a=real(aPer(:,1:N+1));

end




