clear

DO_EXPORT = 0;
nx = 500;
ny = 500;


% xx_b = [0,5];    % xi-coordinate of "begining of" edge
% h = [1,.5,1]; % plateau levels
% % h = [.5,1,.5]; % plateau levels
% theta = [.5,.5]*pi/2; % slope angles (positive values)
% rangeZeta = [-5,10,-pi,1]; % plot range
% exportPath = './SCnumStep2x45deg'; 


% xx_b = [0];    % xi-coordinate of "begining of" edge
% h = [.5,1,]; % plateau levels
% theta = [1]*pi/2; % slope angles (positive values)
% rangeZeta = [-5,10,-pi,1]; % plot range

% xx_b = [0,8,14];    % xi-coordinate of "begining of" edge
% h = [.5,1.5,0.75,1.5]; % plateau levels
% theta = [.5,1.25,1.25]*pi/2; % slope angles (positive values)
% rangeZeta = [-5,18,-pi+.001,1]; % plot range
% exportPath = './SCnumStepMulti'; 
% % rangeZeta = [-5,18,-pi+.1,1]; % plot range
% % exportPath = './SCnumStepMultiSmooth'; 



xx_b = [0,nan];  %2.7726  % xi-coordinate of "begining of" edge
h = [5,2.4,1]; % plateau levels
theta = [1,.5]*pi/2; % slope angles (positive values)
rangeZeta = [-5,50,-pi+.001,1]; % plot range
exportPath = './SCnumBeachRamp'; 
xx_b(2) = xx_b(1) - pi/theta(2)*log(h(3)/h(2));


zzRoots = [xx_b-1i*pi;xx_b+pi./theta.*log(h(2:end)./h(1:end-1))-1i*pi].';


xx = linspace(rangeZeta(1),rangeZeta(2),nx);
yy = linspace(rangeZeta(4),rangeZeta(3),ny)'; % integrating top-to-bottom
% yy = linspace(rangeZeta(3),rangeZeta(4),ny)';  % integrating bottom-to-top
zz = xx + 1i*yy;

% function df
% df = 1;
% for i = 1:length(xx_b)
%     c = (h(i+1)/h(i));%^abs(pi/theta(i)/2);
%     df = df.*h(i)/pi.*((exp(zz-xx_b(i))+1)./(exp(zz-xx_b(i))+c^2)).^(theta(i)/pi);
% end

h = shiftdim(h,-1); theta = shiftdim(theta,-1);xx_b = shiftdim(xx_b,-1);
df = prod(((exp(zz-xx_b)+1)./(exp(zz-xx_b)+(h(2:end)./h(1:end-1)).^(pi/theta))).^(theta/pi),3); %  df ~ 1/tau


df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
z = cumsum([0,df_xh.*diff(xx)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*yy)],1);

% z0 = interp2(xx,yy,real(z),0,0)+1i*interp2(xx,yy,imag(z),0,0);
z0 = interp2(xx,yy,real(z),0,rangeZeta(3))+1i*interp2(xx,yy,imag(z),0,0);
z = z-z0; % orientation constant
K = 1/(-min(imag(z(ny,:))))*max(h);% scaling constant
z = K*z;
df = K*df;

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


% return
hfzz = figure('color','w');
contour(real(zz),imag(zz),real(z),20,'r','linewidth',1);hold on
contour(real(zz),imag(zz),imag(z),10,'b','linewidth',1);
% plot(zetaS,'-ok','linewidth',1)
xlabel('\xi');ylabel('i \sigma');
axis equal tight 
box off
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])

plot(zzRoots,'k*')

if DO_EXPORT
    export_fig(hfz,[exportPath,'_z'],'-pdf','-m2');
    export_fig(hfzz,[exportPath,'_zz'],'-pdf','-m2');
%     print(hfz, '-dpdf','-bestfit',[exportPath,'_z.pdf']); % supports transparancy
%     print(hfzz, '-dpdf','-bestfit',[exportPath,'_zz.pdf']); 
end

