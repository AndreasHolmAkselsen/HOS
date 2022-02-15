clear
g = 9.81;
h = 1;
d = .5;


a = 10+h;
b = 10+h+d;

limZeta = [-5.-b,-a+5,0,2]*h;
limZeta = [a-5,b+5,0,2]*h;

limZeta = [-5,+5,a-1,b+1]*h;



C=0;
K = -2*d/pi;
% dff =@(zz) K*sqrt(zz-ic).*sqrt(zz+ic)./(sqrt(zz-ia).*sqrt(zz+ia).*sqrt(zz-ib).*sqrt(zz+ib));
% dff =@(zz) K*sqrt(zz-b).*sqrt(zz+a)./(zz.*sqrt(zz+b).*sqrt(zz-a));

dff =@(zz) K./sqrt(zz+b).*sqrt(zz+a)./zz.*sqrt(zz-a)./sqrt(zz-b);




nLines = 10;
nPoints = 200;




dzz = 0.01*h; 
zzr = limZeta(1):dzz:limZeta(2);
zzi = flipud((limZeta(3):dzz:limZeta(4)).');
% zzi = (limZeta(3):dzz:limZeta(4)).';
zz = zzr+1i*zzi;
[ny,nx] = size(zz);

df = dff(zz);

df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
f = cumsum([0,df_xh.*diff(zzr)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*zzi)],1) + C; % integrating bottom-to-top

% f = f-mean(f);

% f = 1i*(d-h)+  2*d/pi*( sqrt(zz/d).*sqrt(1+zz/d)-log(sqrt(zz/d)+sqrt(1+zz/d)) ) ;


hf = figure('color','w');
contour(real(f),imag(f),real(zz),20,'b');hold on
contour(real(f),imag(f),imag(zz),20,'r');

return
