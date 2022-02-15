clear
H = 1;
D = .5;

nLines = 10;
nPoints = 300;

k0 = 2*pi/(H);

% limZz = [-7.5,15,0,7]*H;
% f = @(zz) D/pi*( sqrt(zz-1).*sqrt(zz+1)-2*asinh(sqrt(.5*(zz-1))) ) - 1i*(H-D);

limZz = [-3,3,0,1.5]*H;
f = @(zz) 1i*(D-H)+  2*D/pi*( sqrt(zz/D).*sqrt(1+zz/D)-log(sqrt(zz/D)+sqrt(1+zz/D)) ) ;

% limZz = [-3,3,0,1.5]*H;
% f = @(zz) 1i*(D-H) + D/pi*( sqrt(zz/D+1).*sqrt(zz/D-1) - log( zz/D + sqrt(zz/D+1).*sqrt(zz/D-1) )) ;


hfz = figure('color','w'); haz = gca;
% contour(real(f(zz)),imag(f(zz)),real(zz));hold on
% contour(real(f(zz)),imag(f(zz)),imag(zz));
zz = linspace(limZz(1),limZz(2),nLines) + 1i*linspace(limZz(3),limZz(4),nPoints)';
plot(f(zz),'r'); hold on;

zz =  linspace(limZz(1),limZz(2),nPoints) + 1i*linspace(limZz(3),limZz(4),nLines)';
plot(f(zz).','b')
axis equal

plot(real(f(limZz(1:2))),[0,0],'k','linewidth',1)


% map zzS onto the real equidistant points of the line xS
xLR =  real(f(limZz(1:2)));
% xS = linspace(real(f(limZz(1))), real(f(limZz(2))),nLines);
% zz = linspace(limZz(1),limZz(2),100) + 1i*linspace(limZz(3),limZz(4),100)';
% si = scatteredInterpolant(  real(f(zz(:))), imag(f(zz(:))) , zz(:));
% zzS = si(xS, 0*xS );
% plot(real(f(zzS)),imag(f(zzS)),'ok','linewidth',1)
hp = patch([xLR(1),xLR(1),0,0,xLR(end),xLR(end)],[-1.2*H,-H,-H,-H+D,-H+D,-1.2*H],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
xlabel('x');ylabel('i y');
axis tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
box off

hfzz = figure('color','w'); hazz = gca;
zz = linspace(limZz(1),limZz(2),100) + 1i*linspace(limZz(3),limZz(4),100)';
z = f(zz);
contour(real(zz),imag(zz),real(z),'r');hold on
contour(real(zz),imag(zz),imag(z),'b');
% plot(zzS,'-ok','linewidth',1)
axis equal 
box off
xlabel('\xi');ylabel('i \sigma');
set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])


HMax = .3*H;
nArrX = 300;
nArrY = 200;
sigmaUpper = (HMax/2+H)*pi/2;    % from xi->+inf limit
sigmaLower = (-HMax/2+H-D)*pi/2; % from xi->-inf limit

zzArr = linspace(limZz(1),limZz(2),nArrX) + 1i*linspace(sigmaLower,sigmaUpper,nArrY)';
zArr = f(zzArr);
finvIp = scatteredInterpolant(  real(zArr(:)), imag(zArr(:)) , zzArr(:),'linear','none');
plot(hazz,[zzArr(1,1),zzArr(1,end),zzArr(end,end),zzArr(end,1),zzArr(1,1)],'-g')
plot(haz,  [ zArr(1,:),nan, zArr(:,end).',nan,zArr(end,:),nan,zArr(:,1).'],'-g'  )

% zj = linspace(xLR(1),xLR(2),40)  + .5i*HMax*linspace(-1,1,8).';
% zzj = finvIp(real(zj), imag(zj) );
% plot(hazz,zzj,'.k')
% plot(haz,f(zzj),'.k')

xSj = linspace(xLR(1),xLR(2),200);
hj  =  .5*HMax*cos(k0*xSj + 30*pi/180);
zzSj = finvIp(xSj, hj );

plot(haz,xSj,hj,'-k');
plot(hazz,zzSj,'-k');

return
export_fig(hfz,'SC_step','-pdf');
export_fig(hfzz,'SC_step_inv','-pdf');

