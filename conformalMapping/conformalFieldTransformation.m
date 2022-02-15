close all

k = -1;
a = .1;
L = 2*pi;
h = 1;

k0 = pi/2;

A = .1*exp(1i*30*pi/180 );
L = pi;
h = 1;


heta = shiftdim( .75*(A+A') ,-1);
k = shiftdim( [1,-1]*k0  ,-1);

% figure;
% [x,y]= meshgrid(linspace(-L,L,10),linspace(-h,h,100));
% z = x+1i*y; zeta = z + a.*exp(1i*k*z);
% plot(zeta,'r'); hold on;
% 
% [x,y]= meshgrid(linspace(-L,L,100),linspace(-h,h,10));
% z = x+1i*y; zeta = z + a.*exp(1i*k*z);
% plot(zeta.','b')
% axis equal


f =@(zeta) zeta + 1i*sum( heta.*exp(1i*k.*zeta) ,3);
eta = sum( heta.*exp(1i*k.*x) ,3);

x = linspace(-L,L,100);
y = linspace(-h,h,100);
[xx,yy]= meshgrid(x,y);
% z = x+1i*y; zeta = z + a.*exp(1i*k*z);
hf = figure('color','w');
subplot(2,1,1)
contour(real(f(xx+1i*yy)),imag(f(xx+1i*yy)),xx,'r'); hold on;
contour(real(f(xx+1i*yy)),imag(f(xx+1i*yy)),yy,'b');
plot(x, eta,'k','linewidth',2)
axis equal
grid on
ylim([-h,h])
xlabel('x'); ylabel('i y');
%plot inverse
subplot(2,1,2);
contour(xx,yy,real(f(xx+1i*yy)),'r'); hold on;
contour(xx,yy,imag(f(xx+1i*yy)),'b');
contour(xx,yy,imag(f(xx+1i*yy)),[0,0],'k','linewidth',2)
axis equal
grid on
xlabel('\xi'); ylabel('i \sigma');


return
export_fig(hf,'conformal','-pdf')