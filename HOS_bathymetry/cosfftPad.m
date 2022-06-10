function y_p = cosfftPad(y,n_p)
% precise zero-padding/un-padding of high frequencies of sine or cosine fft signal y.
% y: Fourier-space array (y=cosfft(x) or y=sinfft(x)) with fft taken along columns
% n_p: Number of return array spacings (n_p=size(y_p,1)-1)
%
% % example:
% N = 2^4;
% Nd = 3*N;
% ii = (0:N)';
% x = (rand(N+1,1)-.5) +1i*(rand(N+1,1)-.5);
% %%% using cosine
% y = cosfft(x);
% xPad = icosfft(cosfftPad(y,Nd));
% xUnPad = icosfft(cosfftPad(cosfft(xPad),N));
% %%% using sin
% % x([1,end],:) = 0;
% % y = sinfft(x);
% % xPad = isinfft(cosfftPad(y,Nd));
% % xUnPad = isinfft(cosfftPad(sinfft(xPad),N));
% %%%
% figure;
% subplot(211);plot(ii,real(x),'.-',(0:Nd)'*N/Nd,real(xPad),'.--',ii,real(xUnPad),':','linewidth',1); title('real');grid on;
% subplot(212);plot(ii,imag(x),'.-',(0:Nd)'*N/Nd,imag(xPad),'.--',ii,imag(xUnPad),':','linewidth',1); title('imag');grid on;

% Using n = length(y)-1, n_p = length(y_p)-1
n = size(y,1)-1;
m = size(y,2);
if n_p>n
    y_p = n_p/n*[y(1:end-1,:);.5*y(end,:);zeros(n_p-n,m)];
elseif n_p<n
    y_p = n_p/n*[y(1:n_p,:);2*y(n_p+1,:)];
else
    y_p = y;
end

% % Using n = length(y), n_p = length(y_p), n = size(y,1);
% m = size(y,2);
% if n_p>n
%     y_p = (n_p-1)/(n-1)*[y(1:end-1,:);.5*y(end,:);zeros(n_p-n,m)];
% elseif n_p<n
%     y_p = (n_p-1)/(n-1)*[y(1:n_p-1,:);2*y(n_p,:)];
% else
%     y_p = y;
% end

