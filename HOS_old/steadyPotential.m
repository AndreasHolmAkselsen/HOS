clear
close all

x = 0:.01:10;
y = (-3:.01:3)';

% % using two vortices
% zS = [0];% source locations
% Q = [0];% source strengths
% zV = [1-1i,4-4i]; % vortex locations
% Gamma = [1,-3];    % vortex strengths

% using a vortex and a source
zS = [-.5-3i];% source locations
Q = [3];% source strengths
zV = [1-1i]; % vortex locations
Gamma = [1];    % vortex strengths


z = x+1i*y;

zV = shiftdim(zV,-1);
zS = shiftdim(zS,-1);
Q = shiftdim(Q,-1);
Gamma = shiftdim(Gamma,-1);
w = sum(Gamma./(2i*pi).*(log(z-zV) - log(z-conj(zV))),3)...
    + sum(Q./(2*pi).*(log(z-zS) + log(z-conj(zS))),3);

% w = 0*z;
% for i=1:length(zV), w = w + Gamma(i)./(2i*pi).*(log(z-zV(i)) - log(z-conj(zV(i)))); end
% for i=1:length(zS), w = w + Q(i)./(2*pi).*(log(z-zS(i)) + log(z-conj(zS(i)))); end


zq = linspace(x(1),x(end),10) + 1i*linspace(y(1),y(end),10)';
U = conj(  sum(Gamma./(2i*pi).*(1./(zq-zV) - 1./(zq-conj(zV))),3) ...
         + sum(Q/(2*pi).*(1./(zq-zS) + 1./(zq-conj(zS))),3) );


figure('color','w'); 
plot(x([1,end]),[0,0],'k'); hold on;axis equal
contour(real(z),imag(z),imag(w),20);
quiver(real(zq),imag(zq),real(U)./abs(U),imag(U)./abs(U),.25)
