clear

DO_EXPORT = 1;
nx = 3000;
ny = 500;


% xx_b = [0,5];    % xi-coordinate of "begining of" edge
% H = [1,.5,1]; % plateau levels
% % H = [.5,1,.5]; % plateau levels
% theta = [.5,.5]*pi/2; % slope angles (positive values)
% rangeZeta = [-5,10,-pi,1]; % plot range
% exportPath = './SCnumStep2x45deg'; 


% xx_b = [0];    % xi-coordinate of "begining of" edge
% H = [.5,1,]; % plateau levels
% theta = [1]*pi/2; % slope angles (positive values)
% rangeZeta = [-5,10,-pi,1]; % plot range
% 
% xx_b = [0,8,14];    % xi-coordinate of "begining of" edge
% H = [.5,1.5,0.75,1.5]; % plateau levels
% theta = [.5,1.25,1.25]*pi/2; % slope angles (positive values)
% rangeZeta = [-5,18,-pi+.001,1]; % plot range
% exportPath = './SCnumStepMulti'; 
% % rangeZeta = [-5,18,-pi+.1,1]; % plot range
% % exportPath = './SCnumStepMultiSmooth'; 

%%%%%%%%%%  beach tests, 5m depth, 1m deep extension
% % 45 deg
% xx_b = [0,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,2.4,1]; % plateau levels
% theta = [1,.5]*pi/2; % slope angles (positive values)
% rangeZeta = [-5,50,-pi+.001,1]; % plot range
% exportPath = './SCnumBeachRamp'; 

% % 30 deg
% xx_b = [0,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,2.25,1]; % plateau levels
% theta = [90,30]*pi/180; % slope angles (positive values)
% rangeZeta = [-5,50,-pi+.001,1]; % plot range
% exportPath = './SCnumBeachRamp30'; 
%%%%%%%%%%%%

% %%%%%%%%%  beach tests, 5m depth, 2m deep extension sticking out
% % % 45 deg
% % xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% % H = [5,4.7,3.5,2]; % plateau levels
% % theta = [90,180,45]*pi/180; % slope angles (positive values)
% % rangeZeta = [-3.6,10,-pi+.00,1]; % plot range
% % exportPath = './beachRamp_h5D3th45'; 
% 
% % 30 deg
% nx = 3000;
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,4.85,3.33,2]; % plateau levels
% theta = [90,180,30]*pi/180; % slope angles (positive values)
% rangeZeta = [-3.6,15,-pi+.00,1]; % plot range
% exportPath = './beachRamp_h5D3th30'; 
% %%%%%%%%%%%%


% % 30 deg outgoing beach, 10m depth, 1m deep extension  
% nx = 3000;
% D = 3;
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [10,7,3.5,2]; % plateau levels
% theta = [90,180,30]*pi/180; % slope angles (positive values)
% rangeZeta = [-5,14,-pi+.00,1]; % plot range
% exportPath = './SCnumBeachRamp'; 

%%%%%%%%%%%%%%%%%%%% 10m depth, 2m deep extension sticking out
% 30 deg
% nx = 3000;
% D = 4;
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [10,9.1,4.6,2]; % plateau levels
% theta = [90,180,30]*pi/180; % slope angles (positive values)
% rangeZeta = [-5,14,-pi+.00,1]; % plot range
% exportPath = './beachRamp_h10D4th30'; 

% 45 deg
% nx = 3000;
% D = 4;
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [10,8.3,4.9,2]; % plateau levels
% theta = [90,180,45]*pi/180; % slope angles (positive values)
% rangeZeta = [-5,14,-pi+.00,1]; % plot range
% exportPath = './beachRamp_h10D4th45'; 

% % 90 deg
% nx = 3000;
% D = 4;
% xx_b = [0];  %2.7726  % xi-coordinate of "begining of" edge
% H = [10,2]; % plateau levels
% theta = [90]*pi/180; % slope angles (positive values)
% rangeZeta = [-5,14,-pi+.00,1]; % plot range
% exportPath = './beachRamp_h10D4th90'; 

%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%  beach tests, 5m depth, 2m deep extension sticking out
% 30 deg
nx = 3000;
D = 4;
xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
H = [5,4.95,4.33,2]; % plateau levels
theta = [90,180,30]*pi/180; % slope angles (positive values)
rangeZeta = [-3.6,15,-pi+.00,1]; % plot range
exportPath = './beachRamp_h5D4th30'; 

% % 45 deg
% nx = 3000;
% D = 4;
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,4.95,4.43,2]; % plateau levels
% theta = [90,180,45]*pi/180; % slope angles (positive values)
% rangeZeta = [-3.6,10,-pi+.00,1]; % plot range
% exportPath = './beachRamp_h5D4th45'; 

% % 90 deg
% nx = 3000;
% D = 4;
% xx_b = [0];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,2]; % plateau levels
% theta = [90]*pi/180; % slope angles (positive values)
% rangeZeta = [-5,14,-pi+.00,1]; % plot range
% exportPath = './beachRamp_h5D4th90'; 
%%%%%%%%%%%%




for iNan = find(isnan(xx_b))
    xx_b(iNan) = xx_b(iNan-1) - pi./theta(iNan).*log(H(iNan+1)./H(iNan));
end
zzRoots = [xx_b-1i*pi;xx_b+pi./theta.*log(H(2:end)./H(1:end-1))-1i*pi].';


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

H = shiftdim(H,-1); theta = shiftdim(theta,-1);xx_b = shiftdim(xx_b,-1);
df = prod(((exp(zz-xx_b)+1)./(exp(zz-xx_b)+(H(2:end)./H(1:end-1)).^(pi/theta))).^(theta/pi),3); %  df ~ 1/tau


df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
z = cumsum([0,df_xh.*diff(xx)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*yy)],1);

% z0 = interp2(xx,yy,real(z),0,0)+1i*interp2(xx,yy,imag(z),0,0);
z0 = interp2(xx,yy,real(z),0,rangeZeta(3))+1i*interp2(xx,yy,imag(z),0,0);
z = z-z0; % orientation constant
K = 1/(-min(imag(z(ny,:))))*max(H);% scaling constant
z = K*z;
df = K*df;

% plot z-plane
hfz = figure('color','w');hold on
hfz.Position(3) = 750;


% contour(real(z),imag(z),real(zz),20,'r','linewidth',1);
% contour(real(z),imag(z),imag(zz),10,'b','linewidth',1);
plot(z(:,1:round(end/30):end),'r','linewidth',1);  plot(z(1:round(end/10):end,:).' ,'b')
minIz = min(imag(z(ny,:)));
hp = patch(real(z(ny,[1,1:nx,nx])),[1.1*minIz,imag(z(ny,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none'); %,'FaceAlpha',.5


xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
box off


% Lader tank parameters
if exist('D','var')
    D__L = .353;
    L = D/D__L;
    xB = linspace(-L,.25*L,100);
    zB = -D*(xB/L).^2;
    x0 = -interp1(zB(xB<0),xB(xB<0),-2);
    plot(x0+xB,zB,'--k','LineWidth',2)
end

% return
hfzz = figure('color','w');
hfzz.Position(3) = 750;
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

