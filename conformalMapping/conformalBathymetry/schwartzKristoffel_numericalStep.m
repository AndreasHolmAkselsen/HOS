clear
% global h1 h2



h_s = 1;
h_d = .5;
nx = 300;
ny = 300;

limZeta = [-15*h_s,15*h_s,-pi,1];

theta = .5*pi/2;

xx = linspace(limZeta(1),limZeta(2),nx);
yy = linspace(limZeta(4),limZeta(3),ny)'; % integrating top-to-bottom
% yy = linspace(limZeta(3),limZeta(4),ny)';  % integrating bottom-to-top
zz = xx + 1i*yy;


% c = (h_d/h_s).^(pi/theta/2);
% df = h_d/pi.*sqrt((exp(zz)+1)./(exp(zz)+c^2));
% df = h_d/pi.*((exp(zz)+1)./(exp(zz)+c^2)).^(theta/pi);
tau = ((exp(zz)+(h_d/h_s).^(pi/theta))./(exp(zz)+1)).^(theta/pi);
df = h_d./(pi.*tau);

df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
z = cumsum([0,df_xh.*diff(xx)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*yy)],1);

% z0 = interp2(xx,yy,real(z),0,0)+1i*interp2(xx,yy,imag(z),0,0);
z0 = interp2(xx,yy,real(z),0,-pi)+1i*interp2(xx,yy,imag(z),0,0);
z = z-z0;

% plot z-plane
hfz = figure('color','w');hold on

contour(real(z),imag(z),real(zz),20,'r','linewidth',1);
contour(real(z),imag(z),imag(zz),10,'b','linewidth',1);
minIz = min(imag(z(ny,:)));
hp = patch(real(z(ny,[1,1:nx,nx])),[1.1*minIz,imag(z(ny,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none'); %,'FaceAlpha',.5


xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
box off

hfzz = figure('color','w');
contour(real(zz),imag(zz),real(z),20,'r','linewidth',1);hold on
contour(real(zz),imag(zz),imag(z),10,'b','linewidth',1);
% plot(zetaS,'-ok','linewidth',1)
xlabel('\xi');ylabel('i \sigma');
axis equal tight 
box off
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])

return
export_fig(hfz,'SCnumStep45deg_z','-pdf');
export_fig(hfzz,'SCnumStep45deg_zz','-pdf');

% 
% set(hfz,'Units','inches');
% screenposition = get(hfz,'Position');
% set(hfz,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print(hfz, '-dpdf','-painters',['SCnumStep45deg_z.pdf']); % supports transparancy


% print(hfz, '-dpdf','-bestfit',['SCnumStep45deg_z.pdf']); % supports transparancy
% print(hfzz, '-dpdf','-bestfit',['SCnumStep45deg_zz.pdf']);


