rng(1)

% % x = tanh(2*ii/N)-.4*cos(8*ii/N*pi).*exp(ii/N)+1;

%% Using N = length(x)-1, Nd = length(xPad)-1
N = 2^4;
Nd = 3*N;
ii = (0:N)';
x = (rand(N+1,1)-.5) +1i*(rand(N+1,1)-.5);
%%% using cosine
y = cosfft(x);
xPad = icosfft(cosfftPad(y,Nd));
xUnPad = icosfft(cosfftPad(cosfft(xPad),N));
%%% using sin
% x([1,end],:) = 0;
% y = sinfft(x);
% xPad = isinfft(cosfftPad(y,Nd));
% xUnPad = isinfft(cosfftPad(sinfft(xPad),N));
%%%
figure;
subplot(211);plot(ii,real(x),'.-',(0:Nd)'*N/Nd,real(xPad),'.--',ii,real(xUnPad),':','linewidth',1); title('real');grid on;
subplot(212);plot(ii,imag(x),'.-',(0:Nd)'*N/Nd,imag(xPad),'.--',ii,imag(xUnPad),':','linewidth',1); title('imag');grid on;


%% Using N = length(x), Nd = length(xPad)
% N = 2^4;
% Nd = 4*(N-1)+1;
% 
% ii = (0:N-1)';
% x = rand(N,1);
% 
% y = cosfft(x);
% xPad = icosfft(cosfftPad(y,Nd));
% xUnPad = icosfft(cosfftPad(cosfft(xPad),N));
% 
% figure;
% subplot(211);plot(ii,real(x),'.-',(0:Nd-1)'*(N-1)/(Nd-1),real(xPad),'.--',ii,real(xUnPad),':','linewidth',1); title('real');grid on;%
% subplot(212);plot(ii,imag(x),'.-',(0:Nd-1)'*(N-1)/(Nd-1),imag(xPad),'.--',ii,imag(xUnPad),':','linewidth',1); title('imag');grid on;%





% a = sum([.5*A(1);A(2:end-1);.5*A(end)].'.*cos(k.*x),2)./(n-1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

% function f_AA=zeroPaddingCosModes(f,Nd)
% A=getCosModeAmplitudes(f);
% [m,N]=size(A);
% f_AA=getFunctionFromCosModes(Nd/(N-1)*[A,zeros(m,Nd-N+1)]);
% end
% function A=getCosModeAmplitudes(a)
% N=size(a,2)-1;
% aPer=[a,fliplr(a(:,2:end-1))];
% A=fft(aPer,size(aPer,2),2);
% A=real(A(:,1:N+1));
% end
% 
% function a=getFunctionFromCosModes(A)
% N=size(A,2)-1;
% aPer=ifft([A,fliplr(A(:,2:end-1))],2*N,2);
% a=aPer(:,1:N+1);
% end


% N = 2^3+1;
% Nd = 3*N;
% % example signal from gaussian spectrum
% y = ifftshift(exp(-linspace(-2,2,N)'.^2+2i*pi*rand(N,1)));
% x = ifft(y);
% xPad = ifft(fftPad(y,Nd));
% xUnPad = ifft(fftPad(fft(xPad),N));
% figure;
% subplot(211);plot(0:N-1,real(x),'.-',(0:Nd-1)*N/Nd,real(xPad),'.--',0:N-1,real(xUnPad),':','linewidth',1); title('real')
% subplot(212);plot(0:N-1,imag(x),'.-',(0:Nd-1)*N/Nd,imag(xPad),'.--',0:N-1,imag(xUnPad),':','linewidth',1); title('imag')
