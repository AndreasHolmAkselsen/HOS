clear
g = 9.81;
h = 1;
d = .5;
t = 0;

nLines = 10;
nPoints = 300;


C=0;
K = 2*d/pi;

% df =@(zz) K*sqrt(zz)./sqrt(zz+d);

limZeta = [-3,3,0.001,1]*h;

dzz = 0.01*h; 
zzr = limZeta(1):dzz:limZeta(2);
zzi = flipud((limZeta(3):dzz:limZeta(4)).');
zz = zzr+1i*zzi;
[ny,nx] = size(zz);


df = K*sqrt(zz)./sqrt(zz+d);
% df = K*sqrt(zz./(zz+d));

df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
f = cumsum([0,df_xh.*diff(zzr)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*zzi)],1) + C; % integrating bottom-to-top


% f = 1i*(d-h)+  2*d/pi*( sqrt(zz/d).*sqrt(1+zz/d)-log(sqrt(zz/d)+sqrt(1+zz/d)) ) ;


hf = figure('color','w');
contour(real(f),imag(f),real(zz),50,'b');hold on
contour(real(f),imag(f),imag(zz),50,'r');
axis equal
return
