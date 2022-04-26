clear
% global h1 h2



h_d = .3;
h_s = 1.5;
nLines = 10;
nx = 300;%500;
ny = 300;



theta = .5*pi/2;
c = (h_s/h_d).^(pi/theta/2);

K = h_s/pi;
% dff =@(zz) K.*sqrt((exp(zz)+1)./(exp(zz)+c^2));
dff =@(zz) K.*((exp(zz)+1)./(exp(zz)+c^2)).^(theta/pi);
limZeta = [-15*h_d,30*h_d,-pi,1];


xx = linspace(limZeta(1),limZeta(2),nx);
yy = linspace(limZeta(4),limZeta(3),ny)'; % integrating top-to-bottom
% yy = linspace(limZeta(3),limZeta(4),ny)';  % integrating bottom-to-top
zz = xx + 1i*yy;
df = dff(zz);
df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
z = cumsum([0,df_xh.*diff(xx)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*yy)],1);

z0 = interp2(xx,yy,real(z),0,0)+1i*interp2(xx,yy,imag(z),0,0);
z = z-z0;

% plot z-plane
hfz = figure('color','w');hold on

contour(real(z),imag(z),real(zz),nLines,'r','linewidth',1);
contour(real(z),imag(z),imag(zz),nLines,'b','linewidth',1);
minIz = min(imag(z(ny,:)));
hp = patch(real(z(ny,[1,1:nx,nx])),[1.1*minIz,imag(z(ny,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');


xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
box off

hfzz = figure('color','w');
contour(real(zz),imag(zz),real(z),'r','linewidth',1);hold on
contour(real(zz),imag(zz),imag(z),'b','linewidth',1);

% plot(zetaS,'-ok','linewidth',1)
xlabel('\xi');ylabel('i \sigma');
axis equal tight 
box off
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
return
export_fig(hfz,'CCMei_step','-pdf');
export_fig(hfzz,'CCMei_step_inv','-pdf');
export_fig(hfLam,'CCMei_step_lambda','-pdf');
export_fig(hfT,'CCMei_step_t','-pdf');

% function y = sqrt2(x)
%     y = sqrt(x);
%     ii = imag(x)<0;
%     y(ii) = -y(ii);
% end

% function y = log2(x)
%     y = log(x);
%     ii = imag(x)<0;
%     y(ii) = y(ii) + 2i*pi;
% end
% 
% function fz=ffz(zeta)
% global h1 h2
% c = h2/h1;
% % lambda = exp(zeta); % surface at imag(zeta) = pi, bed at imag(zeta) = 0
% lambda = exp(zeta+1i*pi); % surface at imag(zeta) = 0, bed at imag(zeta) = -pi
% t = sqrt((lambda-c^2)./(lambda-1));
% % t = sqrt2(lambda-c^2)./sqrt2(lambda-1);
% fz = -1i*h1 + h2/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
% end



