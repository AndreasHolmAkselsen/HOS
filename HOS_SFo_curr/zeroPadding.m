function f_AA=zeroPadding(f,Nd)
A=getModeAmplitudes(f);
[m,N]=size(A);
f_AA=getFunctionFromModeAmplitudes(Nd/N*[A,zeros(m,Nd-N)]);
end