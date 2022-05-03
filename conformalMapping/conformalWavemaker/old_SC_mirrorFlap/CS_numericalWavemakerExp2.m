clear
% global h1 h2
close all


wbl = 1.6;
wbOveWater = .0;
thetaMax = 40*pi/180;

limZeta = [0,5*wbl,pi];



nx = 500;
ny_minus = 300;
ny_plus = 0;


% integrate from top-top-bottom, right to left:
xx = linspace(limZeta(1),limZeta(2),nx);
% xx = linspace(limZeta(2),limZeta(1),nx);
% yy = linspace(limZeta(4),limZeta(3),ny)'; 


yy = linspace(0,limZeta(3),ny_minus)';
dy = yy(1)-yy(2);
yy = [(ny_plus:-1:1)'*dy;yy ];
ny = ny_minus+ny_plus;
assert(ny==numel(yy));

yy = flipud(yy);

% yy = yy-pi;

%     theta = thetaMax*cos( 2*pi/T*t );
theta = thetaMax;


[z,zz] = fz(xx,yy,theta,wbl,wbOveWater);



% plot z-plane
hfz = figure('color','w');hold on
contour(real(z),imag(z),real(zz),20,'r','linewidth',1);
contour(real(z),imag(z),imag(zz),10,'b','linewidth',1);
% minIz = min(imag(z(ny,:)));
% hp = patch(real(z(ny,[1,1:nx,nx])),[1.1*minIz,imag(z(ny,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none'); %,'FaceAlpha',.5

% lowerCorner1 = z(end,end);
% xWb = [0,0,-(wbl+wbOveWater)*tan(theta)] + real(lowerCorner1);
% yWb = [imag(lowerCorner1), -wbl,wbOveWater];
% plot(xWb,yWb,'k','LineWidth',2.5);

xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
box off

% 
% % plot oposite angle
% theta = -theta;
% [zz,z] = fz(xx,yy,theta,wbl,wbOveWater);
% contour(real(z),imag(z),real(zz),20,'r','linewidth',1);
% contour(real(z),imag(z),imag(zz),10,'b','linewidth',1);
% minIz = min(imag(z(ny,:)));
% hp = patch(real(z(ny,[1,1:nx,nx])),[1.1*minIz,imag(z(ny,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none'); %,'FaceAlpha',.5
% lowerCorner2 = z(end,end);
% xWb = [0,0,-(wbl+wbOveWater)*tan(theta)] + real(lowerCorner2);
% yWb = [imag(lowerCorner2), -wbl,wbOveWater];
% plot(xWb,yWb,'k','LineWidth',2.5);
% 
% fprintf('Horizontal wavemaker shift: %g\n',lowerCorner1-lowerCorner2);


% hfzz = figure('color','w');
% contour(real(zz),imag(zz),real(z),20,'r','linewidth',1);hold on
% contour(real(zz),imag(zz),imag(z),10,'b','linewidth',1);
% % plot(zetaS,'-ok','linewidth',1)
% xlabel('\xi');ylabel('i \sigma');
% axis equal tight 
% box off
% set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])

return
export_fig(hfz,'SCnumStep45deg_z','-pdf');
export_fig(hfzz,'SCnumStep45deg_zz','-pdf');

function [z,zz] = fz(xx,yy,theta,wbl,wbOveWater)

zz = xx + 1i*yy;

thp = theta/pi;
% c = -(wbl+wbOveWater)*tan(theta) + 1i*wbOveWater;
% d = sqrt(pi).*(wbl-1i*c).*exp(-1i*theta)./(gamma(1+thp)*gamma(.5-thp));

d = 1;
k = 1;%exp(1i*.39*pi);
p = 1.5;

lam = zz;
% lam = (zz-p).^2;

% lam = 1- exp(-zz);
% lam = exp(-zz)-1;
% lam = exp(-1./zz);
% lam =sinh(zz);
% d=1i;c=0;

df = k.*(1-d^2./lam.^2).^thp./sqrt(lam.^2-p^2); % k=1


df_xh = .5*(df(1,1:end-1)+df(1,2:end));
z1 = cumsum([0,df_xh.*diff(xx)],2) ;
df_yh = .5*(df(1:end-1,:)+df(2:end,:));
z =  cumsum([z1;df_yh.*diff(1i*yy)],1);

% df_yh = .5*(df(1:end-1,1)+df(2:end,1));
% z1 = cumsum([0;df_yh.*diff(1i*yy)],1);
% z1 = z1-z1(yy==0) + 1i*wbOveWater;
% df_xh = .5*(df(:,1:end-1)+df(:,2:end));
% z =  cumsum([z1,df_xh.*diff(xx)],2) ;
end
