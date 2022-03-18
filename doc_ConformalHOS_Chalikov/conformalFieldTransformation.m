close all

% k = -1;
% a = .1;
% L = 2*pi;
% h = 1;

k0 = pi/2;

A = .15*exp(1i*30*pi/180 );
L = pi;
h = 1;


heta = .5*[A,A'];
kx = [1,-1]*k0;
kx3 =  shiftdim( kx ,-1);
heta3 = shiftdim(heta ,-1);

% f =@(zeta) zeta + 1i*sum( heta3.*exp(1i*kx3.*zeta) ,3);  % non-decaying
% f =@(zeta) zeta + 1i*conj(A).*exp(-1i*k0.*zeta); % decaying


f =@(zeta) zeta + 2i*sum( conj(heta3).*(kx3>0).*exp(-1i*kx3.*zeta) ,3); % decaying


eta =@(x) sum( heta3.*exp(1i*kx3.*x) ,3);

hf = figure('color','w');
[x,y]= meshgrid(linspace(-L,L,10),linspace(-h,h,100));
plot(f(x+1i*y),'r'); hold on;
[x,y]= meshgrid(linspace(-L,L,100),linspace(-h,h,10));
plot(f(x+1i*y).','b');
x = linspace(-L,L,100);
plot(x, eta(x),'k','linewidth',2)
grid on
axis equal
xlabel('x'); ylabel('i y');
axis([-L,L,-h,h])

% x = linspace(-L,L,100);
% y = linspace(-h,h,100);
% [xx,yy]= meshgrid(x,y);
% % z = x+1i*y; zeta = z + a.*exp(1i*k*z);
% hf = figure('color','w');
% subplot(2,1,1)
% contour(real(f(xx+1i*yy)),imag(f(xx+1i*yy)),xx,'r'); hold on;
% contour(real(f(xx+1i*yy)),imag(f(xx+1i*yy)),yy,'b');
% plot(x, eta(x),'k','linewidth',2)
% axis equal
% grid on
% ylim([-h,h])
% xlabel('x'); ylabel('i y');
% %plot inverse
% subplot(2,1,2);
% contour(xx,yy,real(f(xx+1i*yy)),'r'); hold on;
% contour(xx,yy,imag(f(xx+1i*yy)),'b');
% contour(xx,yy,imag(f(xx+1i*yy)),[0,0],'k','linewidth',2)
% axis equal
% grid on
% xlabel('\xi'); ylabel('i \sigma');


return
export_fig(hf,'conformal','-pdf')