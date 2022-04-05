clear
close all
global h1 h2 Wxx Wx domainType


DO_EXPORT = 0;
exportPath = './figures/';


% domainType = 'simple';
% domainType = 'logstip';
domainType = 'double'; 


% h1 = .25;
% h2 = 1;
h2 = .5;
h1 = 1;
nx = 101;
aEta = .075*max(h1,h2);
kk_cut_factor = 1; % low-pass filtering factor: 1.0: no filters.
nWaves = 2;
exportPrefix = '';

% h1 = .25;
% h2 = 1;
% % h2 = .25;
% % h1 = 1;
% nx = 101;
% aEta = .1*max(h1,h2);
% kk_cut_factor = 1; % low-pass filtering factor: 1.0: no filters.
% nWaves = 4;
% exportPrefix = '';


% h1 = .25;
% h2 = 1;% 
% % h2 = .25;
% % h1 = 1;
% nx = 201;
% aEta = .075*max(h1,h2);
% kk_cut_factor = 1; % low-pass filtering factor: 1.0: no filters.
% nWaves = 4;
% exportPrefix = '';

xL = -2.5*max(h1,h2);
xR = 2.5*max(h1,h2);
Wx = 2.0*max(h1,h2);
nContours = 10;

% assume a surface solution:
L = xR-xL;
k0 = nWaves*2*pi/L;
dx = L/nx;
xj = xL + (0:nx-1)*dx;
% hj  =  A*cos(k0*xj + 30*pi/180);
% phiSj  =  .1*sin(k0*xj + 10*pi/180)+.05*sin(2*k0*xj + 80*pi/180);
hj  =  aEta*cos(k0*xj + 30*pi/180);
phiSj  =  .1*sin(k0*xj + 70*pi/180)+.025*sin(2*k0*xj + 80*pi/180);


% Create inverse function through scattered interpolation:
nArrX = 300;
nArrY = 100;
switch domainType
    case 'simple'
        fz = @fzSimple;
        yyUpper = (max(hj)+h1)*pi/2;   % from xi->+inf limit
        yyLower = (min(hj)+h2)*pi/2; % from xi->-inf limit
        yyPlotUpper = yyUpper;
        yyPlotLower = 0;     
    case {'logstip','double'}
        if strcmp(domainType,'logstip')
            fz = @fzStrip;
        else
            fz = @fzDouble;
            Wxx = fzero( @(Wxx)-2*fzStrip( -Wxx/2 )-Wx,0);
        end
        
        xxL0 = fzero(@(xx) real(fz(xx))-Wx,-1); % numerical inverse, upper left corner of z-domain
        assert(isfinite(xxL0))
        dzdzz_min = min(abs(dfzStrip0(linspace(xxL0,-xxL0,21),h1,h2)));
        yyUpper = 1.4*max(hj)./dzdzz_min;
        yyLower = 1.4*min(hj)./dzdzz_min;
%         yyUpper = 1;
%         yyLower = -1;        
        yyPlotUpper = 1;
        yyPlotLower = -pi;   

    otherwise
        error('Domain type''%s'' not recognised.',domainType);
end


dxxExtra = .1*h1;
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
assert( abs(xxj(1)-xx_ip(1))<1e-9*dx && abs(xxj(nx)-xx_ip(nx))<1e-9*dx)
xx_ip([1,nx]) = xxj([1,nx]);
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
hww(abs(kk)>kk_cut_factor*max(kk)) = 0;
fww = @(zz) sum( hww.*exp(1i*(zz+1i*y0).*kk )./cosh(kk*HScale) ,3);


% % assuming flat surface: % 
% % hww = sum( phiS_ip.*exp(-1i*xxj_ip.*kk ) ,2)/nx;
% hww = shiftdim(fft(phiS_ip),-1).*exp(-1i*kk*xxj(1))/nx;
% fww = @(zz) sum( hww.*exp(1i*(zz+1i*y0).*kk )./cosh(kk*HScale) ,3);
% assert( max(abs( real(fww(xxj_ip))-phiS_ip  )) < 1e-9)


% if DO_PLOT
    hf = figure('color','w','position',[436 63 637 889]); 
    zz = linspace(xxL,xxR,100) + 1i*linspace(yyPlotLower,yyPlotUpper,100)';
    z = fz(zz); 
    ww = fww(zz);
    etaIp = interp1(real(zzSj),imag(zzSj),real(zz(1,:)),'linear','extrap');
    ww(imag(zz)>etaIp) = nan;
    ww( real(zz)<xxL | real(zz)>xxR) = nan;
    phiLim = [max(1.2*min(phiSj),min(real(ww(:)))) , min(1.2*max(phiSj),max(real(ww(:))))];
    phiLevels = linspace(phiLim(1),phiLim(2),nContours);
    
    % plot the z-plane
    haz = subplot(3,1,1); hold on;
    title('z-plane');xlabel('x');ylabel('i y');box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
    [~,hcz] = contourf(real(z),imag(z),real(ww),phiLevels,'LineStyle','none');
    
    zPhi = fz(linspace(xxL,xxR,10) + 1i*linspace(yyPlotLower,yyPlotUpper,200)');
    zPsi = fz(linspace(xxL,xxR,200) + 1i*linspace(yyPlotLower,yyPlotUpper,10)');
    plot(zPhi,'r','linewidth',1); hold on; plot(zPsi.' ,'b')
%     axis equal
    switch domainType
        case {'simple','logstip'}
            hp = patch([real(zPsi(1))*[1,1],0,0,real(zPsi(1,end))*[1,1]],[-1.2*max(h1,h2),-h1,-h1,-h2,-h2,-1.2*max(h1,h2)],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
        case 'double'
            hp = patch([real(zPsi(1))*[1,1],-Wx/2*[1,1],Wx/2*[1,1],real(zPsi(1,end))*[1,1]],[-1.1*max(h1,h2),-h2,-h2,-h1,-h1,-h2,-h2,-1.1*max(h1,h2)],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
    end
    xlabel('x');ylabel('i y');
    axis tight 
    plot(haz,  [ zArr(1,:),nan, zArr(:,end).',nan,zArr(end,:),nan,zArr(:,1).'],'--k','linewidth',2 )
%     plot(haz,xj,hj,'-k','linewidth',2);
    surface(xj.*[1;1],hj.*[1;1],0*[xj;xj],[phiSj;phiSj],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    
    
    % plot the zz-plane
    hazz = subplot(3,1,2); hold on    
    title('\zeta-plane'); xlabel('\xi');ylabel('i \sigma'); box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
    
    [~,hczz] = contourf(real(zz),imag(zz),real(ww),phiLevels,'LineStyle','none');
%     hFills = hczz.FacePrims;
%     [hFills.ColorType] = deal('truecoloralpha');
%     for i=1:length(hFills), hFills(i).ColorData(4)=200; end
    clim0 = hazz.CLim;
    contour(real(zz),imag(zz),real(z),'r','linewidth',1);
    contour(real(zz),imag(zz),imag(z),'b','linewidth',1);
    hazz.CLim = clim0;
    
    plot(hazz,[zzArr(1,1),zzArr(1,end),zzArr(end,end),zzArr(end,1),zzArr(1,1)],'--k','linewidth',2 )
%     plot(hazz,zzSj,'-k','linewidth',2);
    surface(real(zzSj).*[1;1],imag(zzSj).*[1;1],0*[xj;xj],[phiSj;phiSj],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);

    
    % evaluate ww at surface to make sure it matches phiS there.
    haPhiS = subplot(3,1,3); hold on
    plot(xj,phiSj,'r',xj,real(fww(zzSj)),'--b','linewidth',1.5)
    legend({'\phi^S(x)','Re\omega[\zeta^S(x)])'},'fontsize',11)
    xlabel('x'); grid on;
    ylim(phiLim)
    

%     if strcmp(domainType,'simple')
%         [haPhiS.YLim,hazz.CLim,haz.CLim] = deal(phiLim);
%         [hczz.LevelList,hcz.LevelList] = deal(linspace(phiLim(1),phiLim(2),10));
%     end
% end


if DO_EXPORT
    if ~isfolder(exportPath), mkdir(exportPath); end
    exportName = sprintf('%s%s_nWaves%d_h1_%g_h2_%g_nx%d_aEta%.2g_kCutF%.2g_L%.2g',exportPrefix,domainType,nWaves,h1,h2, nx,aEta,kk_cut_factor,L); exportName(exportName=='.')='p';
    copyfile('./conformalSurfaceProjection.m',[exportPath,'/script_',exportName,'.m']) 
    savefig(hf,[exportPath,'/',exportName]);
    export_fig(hf,[exportPath,'/',exportName],'-pdf','-png');
end


function z = fzDouble(zz)
global Wxx Wx
    assert(~any(real(zz)==0,'all'), 'It is assumed that no xi values equal zero.')
    fzL = -fzStrip( -(zz+Wxx/2));
    fzR = +fzStrip( +(zz-Wxx/2));
    z = (fzL-Wx/2).*(real(zz)<0) + (fzR+Wx/2).*(real(zz)>0);  
end

function z = fzSimple(zz)
global h1 h2
    assert(h1>h2,'h1>h2 assumed in the ''simple'' configuration.');
    d = h1-h2;
    z = -1i*h2+  2*d/pi*( sqrt(zz/d).*sqrt(1+zz/d)-log(sqrt(zz/d)+sqrt(1+zz/d)) );
end

% function z = fzStrip(zz)
% global h1 h2
%     c = h1/h2;
%     lambda = exp(zz+1i*pi); % surface at imag(zz) = 0, bed at imag(zz) = -pi
%     t = sqrt((lambda-c^2)./(lambda-1));
%     z = -1i*h2 + h1/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
% end

function z=fzStrip(zeta)
global h1 h2
if h1<h2
    z=fzStrip0(zeta,h1,h2);
else
    z=-fzStrip0(-zeta,h2,h1);
end
end

function z = fzStrip0(zz,h1,h2)
    c = h2/h1;
    lambda = -exp(zz); % surface at imag(zz) = 0, bed at imag(zz) = -pi
    t = sqrt((lambda-c^2)./(lambda-1));
    z = -1i*h1 + h2/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
end

function dzdzz = dfzStrip0(zz,h1,h2)
    c2 = (h2/h1)^2; lambda = -exp(zz); 
    t = sqrt((lambda-c2)./(lambda-1));
    dzdzz =  -lambda*h2/pi.* (c2-1).^2./( (c2-t.^2).*(t.^2-1).*(lambda-1).^2.*sqrt((lambda-c2)./(lambda-1)) );
end


function y = log2(x)
    y = log(x);
    ii = imag(x)<0;
    y(ii) = y(ii) + 2i*pi;
end

