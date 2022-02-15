clear
g = 9.81;
h = 1;
d = .5;


ia = 1i*h;
ib = ia + .5i*h;
ic = ib + 2i*d;

C=0;
K = 2*d/pi;
dff =@(zz) K*sqrt(zz-ic).*sqrt(zz+ic)./(sqrt(zz-ia).*sqrt(zz+ia).*sqrt(zz-ib).*sqrt(zz+ib));


xx = linspace(0,1000*h,500);
figure('color','w');
dfxx = dff(xx);
plot(xx,cumsum(real([0,.5*(dfxx(1:end-1)+dfxx(2:end)).*diff(xx)])),'k')
set(gca,'XScale','log')


nLines = 10;
nPoints = 300;




limZeta = [0.,5,-10,10]*h;

dzz = 0.01*h; 
zzr = limZeta(1):dzz:limZeta(2);
% zzi = flipud((limZeta(3):dzz:limZeta(4)).');
zzi = (limZeta(3):dzz:limZeta(4)).';
zz = zzr+1i*zzi;
[ny,nx] = size(zz);

df = dff(zz);

df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
f = cumsum([0,df_xh.*diff(zzr)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*zzi)],1) + C; % integrating bottom-to-top

% f = f-mean(f);

% f = 1i*(d-h)+  2*d/pi*( sqrt(zz/d).*sqrt(1+zz/d)-log(sqrt(zz/d)+sqrt(1+zz/d)) ) ;


hf = figure('color','w');
contour(real(f),imag(f),real(zz),50,'b');hold on
contour(real(f),imag(f),imag(zz),50,'r');

return
