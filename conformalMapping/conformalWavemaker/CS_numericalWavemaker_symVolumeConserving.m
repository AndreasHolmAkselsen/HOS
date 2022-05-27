clear
% global h1 h2


wbl = 1.6;
wbOveWater = .0;
thetaMax = 40*pi/180;

% limZeta = [0,10*wbl,-2*wbl];
% 
% % h_s = 1;
% % h_d = .5;
% nx = 500;
% ny_minus = 500;
% ny_plus = 100;
% % limZeta = [-15*h_s,15*h_s,-pi,1];
% % theta = .5*pi/2;
% 
% % xx = linspace(limZeta(1),limZeta(2),nx);
% % yy = linspace(limZeta(4),limZeta(3),ny)'; % integrating top-to-bottom
% % % yy = linspace(limZeta(3),limZeta(4),ny)';  % integrating bottom-to-top
% 
% % integrate from top-top-bottom, right to left:
% xx = linspace(limZeta(2),limZeta(1),nx);
% % yy = linspace(limZeta(4),limZeta(3),ny)'; 
% 
% 
% yy = linspace(0,limZeta(3),ny_minus)';
% dy = yy(1)-yy(2);
% yy = [(ny_plus:-1:1)'*dy;yy ];
% ny = ny_minus+ny_plus;
% assert(ny==numel(yy));


nx = 500;
ny = 1000;
xx = linspace(10*wbl,0,nx);
yy = linspace(-0*wbl,6*wbl,ny)';
nlx = 30;
nly = 31;


%     theta = thetaMax*cos( 2*pi/T*t );
theta = -thetaMax;


[zz,z] = fz(xx,yy,theta,wbl,wbOveWater);


% plot z-plane
hfz = figure('color','w','position',[1921 414 1280 603]);hold on
hlr = plot(z(:,round(linspace(1,nx,nlx))),'r','linewidth',1);
hlb = plot(z(round(linspace(1,ny,nly)),:).','b','linewidth',1);

plot(z([1,end],end),'--k')

% lowerCorner1 = z(end,end);
% xWb = [0,0,-(wbl+wbOveWater)*tan(theta)] + real(lowerCorner1);
% yWb = [imag(lowerCorner1), -wbl,wbOveWater];
% plot(xWb,yWb,'k','LineWidth',2.5);

xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
box off


return

% plot oposite angle
theta = -theta;
[zz,z] = fz(xx,yy,theta,wbl,wbOveWater);
hlr = plot(z(:,round(linspace(1,nx,nlx))),'r','linewidth',1);
hlb = plot(z(round(linspace(1,ny,nly)),:).','b','linewidth',1);
% lowerCorner2 = z(end,end);
% xWb = [0,0,-(wbl+wbOveWater)*tan(theta)] + real(lowerCorner2);
% yWb = [imag(lowerCorner2), -wbl,wbOveWater];
% plot(xWb,yWb,'k','LineWidth',2.5);



% hfzz = figure('color','w');
% contour(real(zz),imag(zz),real(z),20,'r','linewidth',1);hold on
% contour(real(zz),imag(zz),imag(z),10,'b','linewidth',1);
% % plot(zetaS,'-ok','linewidth',1)
% xlabel('\xi');ylabel('i \sigma');
% axis equal tight 
% box off
% set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])

return
export_fig(hfz,'dafgz','-pdf');
export_fig(hfzz,'dafgzz','-pdf');

function [zz,z] = fz(xx,yy,theta,wbl,wbOveWater)

zz = xx + 1i*yy;

thp = theta/pi;
% c = -(wbl+wbOveWater)*tan(theta) + 1i*wbOveWater;
% d = sqrt(pi).*(wbl-1i*c).*exp(-1i*theta)./(gamma(1+thp)*gamma(.5-thp));
% df = (1+d^2./zz.^2).^thp; 

h1 = .5 *wbl;
h2 = h1+wbl;
h3 = h2+2*wbl;
h4 = h3 + wbl;
zz2 = zz.^2;
% volume preserving
df = ((zz2+h1^2)./(zz2+h4^2).*((zz2+h3^2)./(zz2+h2^2)).^2).^thp;
% sine, but only one mirror
% df = ( (zz2+h1^2).*(zz2+h3^2)./(zz2+h2^2).^2 ).^thp;
% sine, full mirroring
% df = (1+csch(pi*(zz-1i*h2)./(2*h2)).^2.*sin(pi*wbl./(2*h2)).^2).^thp; 
% H = 2*(h4+h1);
% fcc = @(h) cosh(2*pi/H*zz)-cos(2*pi/H*h);
% df = ((fcc(h3)./fcc(h2)).^2.*fcc(h1)./fcc(h4)).^thp;

df_yh = .5*(df(1:end-1,1)+df(2:end,1));
z1 = cumsum([0;df_yh.*diff(1i*yy)],1);
% z1 = z1-z1(yy==0) + 1i*wbOveWater;
df_xh = .5*(df(:,1:end-1)+df(:,2:end));
z =  cumsum([z1,df_xh.*diff(xx)],2) ;
end
