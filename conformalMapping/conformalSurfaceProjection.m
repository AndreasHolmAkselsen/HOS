clear
close all
global D H Wxx Wx domainType


DO_PLOT = 1;
PLOT_phiS = 1;

% domainType = 'simple';
domainType = 'logstip';
% domainType = 'double';

g = 9.81;
H = 1;
D = .4;

nx = 101;

% xL = -2*H;
% xR = 3*H;




xL = -2.5*H;
xR = 2.5*H;
Wx = 2.5*H;

% assume a surface solution:
L = xR-xL;
k0 = 4*2*pi/L;
dx = L/nx;
xj = xL + (0:nx-1)*dx;
hj  =  .05*H*cos(k0*xj + 30*pi/180);
phiSj  =  .1*sin(k0*xj + 10*pi/180)+.05*sin(2*k0*xj + 80*pi/180);


% Create inverse function through scattered interpolation:
nArrX = 300;
nArrY = 200;
switch domainType
    case 'simple'
        yyUpper = (max(hj)+H)*pi/2;   % from xi->+inf limit
        yyLower = (min(hj)+H-D)*pi/2; % from xi->-inf limit
        yyPlotUpper = yyUpper;
        yyPlotLower = 0;     
        fz = @fzSimple;
    case {'logstip','double'}
        yyUpper = 1;
        yyLower = -1;        
        yyPlotUpper = 1;
        yyPlotLower = -pi;   
        if strcmp(domainType,'logstip')
            fz = @fzStrip;
        else
            fz = @fzDouble;
            Wxx = fzero( @(Wxx)-2*fzStrip( -Wxx/2 )-Wx,0);
        end
    otherwise
        error('Domain type''%s'' not recognised.',domainType);
end


dxxExtra = .1*H;
xxL = fzero(@(xx) real(fz(xx+1i*yyUpper))-xL,xL )-dxxExtra; % numerical inverse, upper left corner of z-domain
xxR = fzero(@(xx) real(fz(xx+1i*yyLower))-xR,xR )+dxxExtra; % numerical inverse, lower left corner of z-domain
zzArr = linspace(xxL,xxR,nArrX) + 1i*linspace(yyLower,yyUpper,nArrY)';
zArr = fz(zzArr);
finvIp = scatteredInterpolant(  real(zArr(:)), imag(zArr(:)) , zzArr(:),'linear','none');

% Find the inverted surface zzSj
zzSj = finvIp(xj, hj );
xxj = real(zzSj).';
yyj = imag(zzSj).';



dxx = diff(real(zzSj([1,end])))./(nx-1);
Lxx = dxx*nx;
dkk = 2*pi/Lxx;
assert(mod(nx,2)==1,'Odd number of points assumed.')
kkp = (1:(nx-1)/2)*dkk;
% kk = shiftdim([0,kkp,-kkp],-1);
kk = shiftdim([0,kkp,-fliplr(kkp)],-1);

xx_ip = xxj(1) + (0:nx-1).'*dxx;
phiS_ip = interp1(xxj,phiSj,xx_ip);
yy_ip = interp1(xxj,yyj,xx_ip);

% map phiS onto the coefficients of ww
switch domainType
    case 'simple'
        HScale = imag(finvIp(0,0));
        y0 = 0;
    case {'logstip','double'}
        HScale = pi;
        y0 = pi;
end
% A = [ones(nx,1),2*cos(kkp.*xxj).*cosh(kkp.*(yyj+y0))./cosh(kkp*HScale),-2*sin(kkp.*xxj).*cosh(kkp.*(yyj+y0))./cosh(kkp*HScale)];
% wReIm = A\(phiSj.');
A = [ones(nx,1),2*cos(kkp.*xx_ip).*cosh(kkp.*(yy_ip+y0))./cosh(kkp*HScale),-2*sin(kkp.*xx_ip).*cosh(kkp.*(yy_ip+y0))./cosh(kkp*HScale)];
tic
wReIm = A\phiS_ip;
fprintf('System solver CPU time: %g\n',toc);
hww = zeros(1,1,nx);
hww(1) = wReIm(1);
hww(2:(nx+1)/2) = wReIm(2:(nx+1)/2) + 1i* wReIm((nx+3)/2:nx);
% hww((nx+3)/2:nx) = conj(hww(2:(nx+1)/2));
hww((nx+3)/2:nx) = conj(hww((nx+1)/2:-1:2));
fww = @(zz) sum( hww.*exp(1i*(zz+1i*y0).*kk )./cosh(kk*HScale) ,3);



% % assuming flat surface:
% dxx = diff(real(zzSj([1,end])))./(nx-1);
% xxj_ip = xxj(1) + (0:nx-1)*dxx;
% phiS_ip = interp1(xxj,phiSj,xxj_ip);
% 
% dkk = 2*pi/(nx*dxx);
% kk  = shiftdim([0:(nx-1)/2, -(nx-1)/2:-1].*dkk,-1);
% 
% HScale = pi;
% y0 = pi;
% % hww = sum( phiS_ip.*exp(-1i*xxj_ip.*kk ) ,2)/nx;
% hww = shiftdim(fft(phiS_ip),-1).*exp(-1i*kk*xxj(1))/nx;
% fww = @(zz) sum( hww.*exp(1i*(zz+1i*y0).*kk )./cosh(kk*HScale) ,3);
% assert( max(abs( real(fww(xxj_ip))-phiS_ip  )) < 1e-9)


% evaluate ww at surface to make sure it matches phiS there.
if PLOT_phiS
    figure('color','w'); plot(xj,phiSj,'r',xj,real(fww(zzSj)),'--b','linewidth',1.5)
end

if DO_PLOT
    hfz = figure('color','w','position',[447 401 637 468]); 
    
    % plot the z-plane
    haz = subplot(2,1,1);
    zPhi = fz(linspace(xxL,xxR,10) + 1i*linspace(yyPlotLower,yyPlotUpper,200)');
    zPsi = fz(linspace(xxL,xxR,200) + 1i*linspace(yyPlotLower,yyPlotUpper,10)');
    plot(zPhi,'r'); hold on; plot(zPsi.' ,'b')
%     axis equal
    switch domainType
        case {'simple','logstip'}
            hp = patch([real(zPsi(1))*[1,1],0,0,real(zPsi(1,end))*[1,1]],[-1.2*H,-H,-H,-H+D,-H+D,-1.2*H],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
        case 'double'
            hp = patch([real(zPsi(1))*[1,1],-Wx/2*[1,1],Wx/2*[1,1],real(zPsi(1,end))*[1,1]],[-1.1*H,-H,-H,-(H-D)*[1,1],-H,-H,-1.1*H],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
    end
    xlabel('x');ylabel('i y');
    axis tight 
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    box off    
    plot(haz,  [ zArr(1,:),nan, zArr(:,end).',nan,zArr(end,:),nan,zArr(:,1).'],'-g','linewidth',1  )
    plot(haz,xj,hj,'-k','linewidth',2);
    title('z-plane')
    
    
    % plot the zz-plane
    hazz = subplot(2,1,2);
    zz = linspace(xxL,xxR,100) + 1i*linspace(yyPlotLower,yyPlotUpper,100)';
    z = fz(zz);   
    contour(real(zz),imag(zz),real(z),':r');hold on
    contour(real(zz),imag(zz),imag(z),':b');
%     axis equal
    box off
    xlabel('\xi');ylabel('i \sigma');
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    
    plot(hazz,[zzArr(1,1),zzArr(1,end),zzArr(end,end),zzArr(end,1),zzArr(1,1)],'-g','linewidth',1 )
    plot(hazz,zzSj,'-k','linewidth',2);
    
    
    ww = fww(zz);
    etaIp = interp1(real(zzSj),imag(zzSj),real(zz(1,:)),'linear','extrap');
    ww(imag(zz)>etaIp) = nan;
    ww( real(zz)<xxL | real(zz)>xxR) = nan;
    contour(real(zz),imag(zz),real(ww),'r');hold on
    contour(real(zz),imag(zz),imag(ww),'b');
    title('\zeta-plane')
end







hfC=figure('color','w');
% switch domainType
%     case 'simple'
%         ii = 1:size(zz,1);
%         jj = 20:size(zz,2);
%     case 'logstip'
%         ii = 1:size(zz,1);
%         jj = 70:size(zz,2);
%     case 'double'   
%         ii = 1:80;
%         jj = 1:size(zz,2);
% end
% contourf(real(zz(ii,jj)),imag(zz(ii,jj)),real(ww(ii,jj))); hold on
contourf(real(zz),imag(zz),real(ww)); hold on
surface(real(zzSj).*[1;1],imag(zzSj).*[1;1],0*[xj;xj],[phiSj;phiSj],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2);
colorbar;
axis equal
box off
xlabel('\xi');ylabel('i \sigma');
set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
title('\zeta-plane')



function z = fzDouble(zz)
global Wxx Wx
    assert(~any(real(zz)==0,'all'), 'It is assumed that no xi values equal zero.')
    fzL = -fzStrip( -(zz+Wxx/2) );
    fzR = +fzStrip( +(zz-Wxx/2) );
    z = (fzL-Wx/2).*(real(zz)<0) + (fzR+Wx/2).*(real(zz)>0);  
end

function z = fzSimple(zz)
global D H
    z = 1i*(D-H)+  2*D/pi*( sqrt(zz/D).*sqrt(1+zz/D)-log(sqrt(zz/D)+sqrt(1+zz/D)) );
end

function z = fzStrip(zz)
global D H
    h_s = H-D;
    c = H/h_s;
    lambda = exp(zz+1i*pi); % surface at imag(zz) = 0, bed at imag(zz) = -pi
    t = sqrt((lambda-c^2)./(lambda-1));
    z = -1i*h_s + H/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
end

function y = log2(x)
    y = log(x);
    ii = imag(x)<0;
    y(ii) = y(ii) + 2i*pi;
end

