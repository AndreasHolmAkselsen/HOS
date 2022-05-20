% close
clear
% global h1 h2


wbl = 1.;
h = 2;
wbOveWater = .25;


dt = 0.1;
tMax = 0;
T = 1;
thetaMax = deg2rad(-4.2);

limZeta = [0,2.5*h,-h];

nx = 300;
ny_minus = 300;
ny_plus = 0;

% integrate from top-top-bottom, right to left:
xx = linspace(limZeta(2),limZeta(1),nx);
% yy = linspace(limZeta(4),limZeta(3),ny)'; 


yy = linspace(0,limZeta(3)-wbOveWater,ny_minus)';
dy = yy(1)-yy(2);
yy = [(ny_plus:-1:1)'*dy;yy ];
ny = ny_minus+ny_plus;
assert(ny==numel(yy));


zz = xx + 1i*yy;


% plot z-plane
hfz = figure('color','w');hold on
xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
box off
% axis([-limZeta(2)-1.2*(wbl+wbOveWater)*tan(thetaMax),0,1.2*yy(end),1.5*yy(1)+wbOveWater])

%     minIz = min(imag(z(ny,:)));
%     hp = patch(real(z(ny,[1,1:nx,nx])),[1.1*minIz,imag(z(ny,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none'); %,'FaceAlpha',.5

hp_floor = plot([0,limZeta(2)],-h*[1,1],'k','LineWidth',2);


% yWb = [imag(z(end)), -wbl,wbOveWater,wbl,-imag(z(end))];
if ny_plus>0
    hp_BM = plot(nan(1,5), [limZeta(3), -wbl,wbOveWater,wbl+2*wbOveWater,-limZeta(3)],'o-k','LineWidth',2.5);
else
    hp_BM = plot(nan(1,3),[limZeta(3), -wbl,wbOveWater],'o-k','LineWidth',2.5);    
end

t = 0;
hp_c=[];
while t<=tMax
    theta = thetaMax*cos( 2*pi/T*t );
    z = fz(xx,yy,theta,h,wbl,wbOveWater);
    
    delete(hp_c)
    [~,hp_c(1)] = contour(real(z),imag(z),real(zz),20,'r','linewidth',1);
    [~,hp_c(2)] = contour(real(z),imag(z),imag(zz),20,'b','linewidth',1);
    
    hp_BM.XData = [0,0,-(wbl+wbOveWater)*tan(theta),zeros(1,2*(ny_plus>0))];
    
    t = t+dt;
    drawnow
end
    
% 
% 
% % plot oposite angle
% theta = -thetaMax;
% z = fz(xx,yy,theta,h,wbl,wbOveWater);
% contour(real(z),imag(z),real(zz),20,'r--','linewidth',1);
% contour(real(z),imag(z),imag(zz),20,'b--','linewidth',1);  
% plot([0,-(wbl+wbOveWater)*tan(theta),0],[-wbl,wbOveWater,wbl+2*wbOveWater],'ok--','LineWidth',2.5);





return
export_fig(hfz,'mapFlap3_positive','-png','-pdf')
export_fig(hfz,'mapFlap3_negative','-png','-pdf')
    
function [z,zz] = fz(xx,yy,theta,h,wbl,wbOveWater)
thp = theta/pi;
% d = sqrt(pi).*(wbl+wbOveWater).*sec(theta)./(gamma(1+thp)*gamma(.5-thp));

zz = xx + 1i*yy;
% df = (1+d^2./zz.^2).^thp; % k=1



H = h + wbOveWater;
D = wbl + wbOveWater;
% d = fixHingeDepth(D/H,theta)*H;


d = fzero(@(d__h) fixHingeDepth(d__h,D/H,theta),D/H)*H;
% fixHingeDepth(d/H,D/H,theta)

% h = -yy(end);
df = (1+csch(pi*zz./(2*H)).^2.*sin(pi*d./(2*H)).^2).^thp; % k=1

df_yh = .5*(df(1:end-1,1)+df(2:end,1));
z1 = cumsum([0;df_yh.*diff(1i*yy)],1);
z1 = z1-z1(yy==0) + 1i*wbOveWater;
df_xh = .5*(df(:,1:end-1)+df(:,2:end));
z =  cumsum([z1,df_xh.*diff(xx)],2) ;

z = z-real(z(end));
end



function err = fixHingeDepth(d__h,wbl__h,theta)

nYInt = 500;

% there's a singularity at zz=-1i*d if theta < 0
% either stop a yy=-d:
% delta_singularity = 1e-3*(theta<0);
%  ... L = sum(.5*(dL(1:end-1)+dL(2:end)) .* diff(yi) ) + delta_singularity ;
% or shift integration path d_xi to the right (into the domain)
d_xi = 1e-3*(theta<0);

yi = linspace(-1,-d__h,nYInt)';
dL = real( (1-csc(pi/2*(yi-1i*d_xi)).^2.*sin(pi/2*d__h).^2).^(theta/pi) );
L = sum(.5*(dL(1:end-1)+dL(2:end)) .* diff(yi) );
err = wbl__h-(1-L);
end
