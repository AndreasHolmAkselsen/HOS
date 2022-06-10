function y_p = fftPad(y,n_p)
% precise zero-padding/un-padding of high frequencies of fft signal y.
% y: Fourier-space array (y=fft(x)) with FFT along columns
% n_p: return array size (n_p=size(y_p,1))
%
% % example:
% N = 2^3+1;
% Nd = 3*N;
% % example signal from gaussian spectrum
% y = ifftshift(exp(-linspace(-2,2,N)'.^2+2i*pi*rand(N,1)));
% x = ifft(y);
% xPad = ifft(fftPad(y,Nd));
% xUnPad = ifft(fftPad(fft(xPad),N));
% figure;
% subplot(211);plot(0:N-1,real(x),'.-',(0:Nd-1)*N/Nd,real(xPad),'.--',0:N-1,real(xUnPad),':'); title('real')
% subplot(212);plot(0:N-1,imag(x),'.-',(0:Nd-1)*N/Nd,imag(xPad),'.--',0:N-1,imag(xUnPad),':'); title('imag')


[n0,m] = size(y);

% case treatment including the -n/2 mode:
if n0<n_p
    if mod(n0,2)==0
        y_p = [y(1:n0/2,:);.5*y(n0/2+1,:);zeros(n_p-n0-1,m);.5*y(n0/2+1,:);y(n0/2+2:n0,:)]*n_p/n0;
%         y_p = [y(1:n0/2,:); 0             ;zeros(n_p-n0-1,m);0            ;y(n0/2+2:n0,:)]*n_p/n0;
    else
        y_p = [y(1:(n0+1)/2,:);zeros(n_p-n0,m);y((n0+1)/2+1:n0,:)]*n_p/n0;
    end
elseif n0>n_p
    if mod(n_p,2)==0
        y_p = [y(1:n_p/2,:);2*real(y(n0-n_p/2+1,:));y(n0-n_p/2+2:n0,:)]*n_p/n0; 
%         y_p = [y(1:n_p/2,:);0               ;y(n0-n_p/2+2:n0,:)]*n_p/n0; 
    else
        y_p = [y(1:(n_p+1)/2,:);y(n0-(n_p+1)/2+2:n0,:)]*n_p/n0; 
    end    
else
    y_p = y;
end


% single-line method that excludes the -n/2 mode:
% assert(mod(n0,2)==0);
% nr = min(n0,n_p);
% y_p = [y(1:nr/2,:);zeros(n_p-nr+1,1);y(n0-nr/2+2:n0,:)]*n_p/n0;  % ignoring highest mode