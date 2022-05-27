clear
% global h1 h2


wbl = 1.6;
h = 1.5*wbl;
wbOverWater = .3;
thetaMax = 30*pi/180;



nx = 500;
ny = 1000;
xx = linspace(10*wbl,0,nx);
yy = linspace(-0*wbl,6*wbl,ny)';
nlx = 30;
nly = 31;

theta = -thetaMax;


[zz,z] = fz(xx,yy,theta,h,wbl,wbOverWater);


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


% return

% plot oposite angle
theta = -theta;
[zz,z] = fz(xx,yy,theta,h,wbl,wbOverWater);
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


function [zz,z] = fz(xx,yy,theta,h,wbl,wbOverWater)

    d = wbl+wbOverWater;
    H = h+wbl+2*wbOverWater;
    [zz,z] = fz0(xx,yy,theta,H,d);
    z = z - 1i*h;
    
end


function [zz,z] = fz0(xx,yy,theta,H,d)

zz = xx + 1i*yy;

% d = wbl;
% H = h+wbl;
sig = H + [-2;-1;1;2]*d;
% sig = h + [-1;0;2;3]*wbl;


% options.Display = 'off';
% sig = newtonraphson(@(sig) fixHingeDepth(sig,theta,H,d),sig,options);
% sig = fminsearch(@(sig) sum(fixHingeDepth(sig,theta,H,d).^2),sig);
sig = LMFnlsq(@(sig) fixHingeDepth(sig,theta,H,d),sig);

zz2 = zz.^2;
df = ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);

df_yh = .5*(df(1:end-1,1)+df(2:end,1));
z1 = cumsum([0;df_yh.*diff(1i*yy)],1);
% z1 = z1-z1(yy==0) + 1i*wbOveWater;
df_xh = .5*(df(:,1:end-1)+df(:,2:end));
z =  cumsum([z1,df_xh.*diff(xx)],2) ;
end



function err = fixHingeDepth(sig,theta,H,d)
d_xi = 1e-6;

df_zz2 = @(zz2) ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);
dy = @(yy) real( df_zz2((1i*yy+d_xi).^2) ); % real beacuse it is integrated by dzz = 1i*dyy;

y = zeros(5,1); sig0 = [0;sig];
for i = 1:4
    y(i+1) = y(i) + integral( dy,sig0(i),sig0(i+1));
end


% err = y(2:5) - (h+[-1;0;2;3]*wbl);
err = y(2:5) - (H+[-2;-1;1;2]*d);

end