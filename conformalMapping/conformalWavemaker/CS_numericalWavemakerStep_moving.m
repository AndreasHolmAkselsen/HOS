close
clear
% global h1 h2


wbl = 1.6;
wbOveWater = .5;


dt = 0.1;
tMax = 0;%10;
T = 1;
thetaMax = deg2rad(20);

limZeta = [0,7*wbl,-2*wbl];

nx = 500;
ny_minus = 500;
ny_plus = 100;

% xx = linspace(limZeta(1),limZeta(2),nx);
% yy = linspace(limZeta(4),limZeta(3),ny)'; % integrating top-to-bottom
% % yy = linspace(limZeta(3),limZeta(4),ny)';  % integrating bottom-to-top

% integrate from top-top-bottom, right to left:
xx = linspace(limZeta(2),limZeta(1),nx);
% yy = linspace(limZeta(4),limZeta(3),ny)'; 


yy = linspace(0,limZeta(3),ny_minus)';
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
axis([-limZeta(2)-1.2*(wbl+wbOveWater)*tan(thetaMax),0,1.2*yy(end),1.5*yy(1)+wbOveWater])

%     minIz = min(imag(z(ny,:)));
%     hp = patch(real(z(ny,[1,1:nx,nx])),[1.1*minIz,imag(z(ny,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none'); %,'FaceAlpha',.5

hp_floor = plot(nan,nan,'k','LineWidth',2);


% yWb = [imag(z(end)), -wbl,wbOveWater,wbl,-imag(z(end))];
yWb = [limZeta(3), -wbl,wbOveWater,wbl+2*wbOveWater,-limZeta(3)];
hp_BM = plot(nan(5,1),yWb,'k','LineWidth',2.5);

theta = thetaMax;
t = 0;
hp_c=[];
while t<=tMax
    theta = thetaMax*cos( 2*pi/T*t );
    z = fz(xx,yy,theta,wbl,wbOveWater);
    
    delete(hp_c)
    [~,hp_c(1)] = contour(real(z),imag(z),real(zz),20,'r','linewidth',1);
    [~,hp_c(2)] = contour(real(z),imag(z),imag(zz),10,'b','linewidth',1);
    
    hp_BM.XData = [0,0,-(wbl+wbOveWater)*tan(theta),0,0] + real(z(end));
    
    hp_floor.XData = real(z(ny,:));
    hp_floor.YData = imag(z(ny,:));
    t = t+dt;
    drawnow
end
    


% plot oposite angle
theta = -theta;
z = fz(xx,yy,theta,wbl,wbOveWater);
contour(real(z),imag(z),real(zz),20,'r--','linewidth',1);
contour(real(z),imag(z),imag(zz),10,'b--','linewidth',1);  
plot([0,-(wbl+wbOveWater)*tan(theta),0] + real(z(end)),[-wbl,wbOveWater,wbl+2*wbOveWater],'k--','LineWidth',2.5);
plot(z(ny,:),'k--','LineWidth',2);



% dtheta = -thetaMax*sin( 2*pi/T*t )* 2*pi/T;
% z_t = fz_t(xx,yy,theta,dtheta,wbl,wbOveWater);
% figure;
% subplot(121);contourf(real(z),imag(z),real(z_t));colorbar;
% subplot(122);contourf(real(z),imag(z),imag(z_t));colorbar;

% % test numerically.
% dt = .02*T;
% z1 = fz(xx,yy,thetaMax*cos( 2*pi/T*t ),wbl,wbOveWater);
% z2 = fz(xx,yy,thetaMax*cos( 2*pi/T*(t+dt) ),wbl,wbOveWater);
% d_t_num = (z2-z1)/dt;
% figure;
% subplot(121);contourf(real(z),imag(z),real(d_t_num));colorbar;
% subplot(122);contourf(real(z),imag(z),imag(d_t_num));colorbar;
% figure;
% subplot(121);contourf(real(z),imag(z),real(z_t-d_t_num),30);colorbar;
% subplot(122);contourf(real(z),imag(z),imag(z_t-d_t_num),30);colorbar;

return
export_fig(hfz,'mapFlap2Th','-png','-pdf')
    
function [z,zz] = fz(xx,yy,theta,wbl,wbOveWater)
thp = theta/pi;
c = -(wbl+wbOveWater)*tan(theta) + 1i*wbOveWater;
d = sqrt(pi).*(wbl-1i*c).*exp(-1i*theta)./(gamma(1+thp)*gamma(.5-thp));

zz = xx + 1i*yy;
df = (1+d^2./zz.^2).^thp; % k=1


df_yh = .5*(df(1:end-1,1)+df(2:end,1));
z1 = cumsum([0;df_yh.*diff(1i*yy)],1);
z1 = z1-z1(yy==0) + 1i*wbOveWater;
df_xh = .5*(df(:,1:end-1)+df(:,2:end));
z =  cumsum([z1,df_xh.*diff(xx)],2) ;
end


function [z_t,zz] = fz_t(xx,yy,theta,dtheta,wbl,wbOveWater)
thp = theta/pi;
c = -(wbl+wbOveWater)*tan(theta) + 1i*wbOveWater;
d = sqrt(pi).*(wbl-1i*c).*exp(-1i*theta)./(gamma(1+thp)*gamma(.5-thp));
zz = xx + 1i*yy;
df_t = dtheta/pi*(1+d^2./zz.^2).^thp.*log(1+d^2./zz.^2);
df_yh = .5*(df_t(1:end-1,1)+df_t(2:end,1));
z1 = cumsum([0;df_yh.*diff(1i*yy)],1);
z1 = z1-z1(yy==0) + 1i*wbOveWater;
df_xh = .5*(df_t(:,1:end-1)+df_t(:,2:end));
z_t =  cumsum([z1,df_xh.*diff(xx)],2) ;
end