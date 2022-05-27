clear 
global timeReached 
timeReached = 0; 

g = 9.81;
% g = 1.0;

%% input
h = 5; % Depth
wbl = 3; % hinge depth
wbOverWater = .5; % flap extention obove quiescent waterline
thetaMax = deg2rad(10);
N_wavelenghtsInL = 20;

% interpolation grid for conformal map
nArrX_far = 200;   % # horizontal points (far field)
nArrX_near = 200;  % # horizontal points (near field)
nArrYDown = 200;    % # vertical points
nTheta = 101;       % # flap angles (odd)

% fz = @fz_yySymmetric;
fz = @fz_volumeConverving;

% for wave breaking example
param.DO_PADDING = 0;
RK4dt = 0;%2.5e-3; % set to zero to use ODE45
%  and  NT_dt =  5 /9/T; lambda = 2*pi; g=1;


relTolODE = 1e-4;% 1e-8;
N_SSGW = 2^12; % number of modes in SSGW solution

% Plot & export options
DO_EXPORT = 0;
EXPORT_MAT = 1;
PLOT_MAP = 0;
exportPrefix = 'BM_';
exportPath = './figures/';
exportFormatsMap = {'-pdf','-png'};
exportFormats = {'-png','-pdf','-m2'};


%% Simulation
nx = 2^10;

% % specifying period
T = 2.0; omega=2*pi/T; 
kTemp = findWaveNumbers(omega,h(1),0,0);
L = 2*pi/kTemp*N_wavelenghtsInL;


flapPlotTimes = linspace(0,T,50); % leave empty skips flap plotting

% Simulation/plotting time
NT_dt = 1;
dt = NT_dt*T;
param.t_end = 9*dt;
    

% TRamp = 0*T;
% param.nonLinRamp = @(t) max(0,1-exp(-(t/TRamp)^2));
% param.iModeCut = 2*(NWaves + (param.M+5));


param.M = 5; 
param.iModeCut = inf;
param.kd__kmax = .5;
param.rDamping = .25;


dx = L/nx;
% x = (0:nx-1)'*dx;
x = linspace(0,L,nx+1)';x(end)=[];
fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-param.iModeCut/ (nx/2),0) )

t0 = 0;
initialStepODE = 1e-3*T;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);%,'MaxStep',1e-4);

fileName = sprintf('%sT%.2f_M%d_h%.2f_wbl%.2f_thetaMax%.1f_L%.3g_dt%.3gT_nx%d_pad%d_ikCut%.4g_Md%.2g_r%.2g',exportPrefix,T,param.M,h,wbl,thetaMax*180/pi,L,NT_dt,nx,param.DO_PADDING,param.iModeCut,param.kd__kmax,param.rDamping); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./proto_beachRamp.m',[exportPath,'/m/',fileName,'.m']) 
end




%%  map preparation
xxIP_near = linspace( 0, 5*wbl, nArrX_near);
xxIP_far = linspace(xxIP_near(end),1.1*L,nArrX_far+1); xxIP_far(1) = []; 
xxIP = [xxIP_near,xxIP_far];

yyUpper = 0;%wbOverWater;
assert(mod(nTheta,2)==1,'Odd number of angles assumed.')
thetaIP = shiftdim(  thetaMax*linspace(-1,1,nTheta),-1);


yyIP = linspace(-h,0,nArrYDown)'; % ensure that we capture the line yy=0
dy = yyIP(2)-yyIP(1);
yyIP = [yyIP;(dy:dy:wbOverWater)'];
assert(yyIP(nArrYDown)==0)

[zIP,f_zIP,zzIP] = fz(xxIP,yyIP,thetaIP,h,wbl,wbOverWater);

% % plot interpolation basis to inspect resolution:
% iThetaToPlot = 1;
% if ~isempty(iThetaToPlot)
%     iPlotEnd = round(size(zIP,2)/10);
%     figure('color','w');hold on; grid on
%     plot(0,-wbl,'ok','markersize',8,'linewidth',2);
% %     axis([-(wbl+wbOverWater)*sin(thetaMax),real(zIP(1,iPlotEnd,1)),-h,wbOverWater])
%     for iPlot = iThetaToPlot
%         if iPlot>1, delete([hlr;hlb]); end
%         hlr = plot(zIP(:,1:round(iPlotEnd/30):iPlotEnd,iPlot),'r','linewidth',1);  
%         hlb = plot(zIP(1:round(size(zIP,1)/10):end,1:iPlotEnd,iPlot).','b','linewidth',1);
%         drawnow
% %         minIz = min(imag(zIP(1,iPlotStart:end,iPlot))); patch(real(zIP(1,[1,iPlotStart:end,iPlotEnd],iPlot)),[1.1*minIz,imag(zIP(1,1:iPlotEnd,iPlot)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
%     end
% end
% % assert(all(real(zIP(nArrYDown,1,:))<x(1)&min(real(zIP(:,end,:)),[],1)>x(end)),'Physical domain [%.3g,%.3g] out of range of interpolation range [%.3g,%.3g]. Extend interpolation range.',x(1),x(end),real(zIP(nArrYDown,1)),real(zIP(nArrYDown,end)))
% assert(all(min(real(zIP(:,end,:)),[],1)>x(end)),'Physical domain out of interpolation range.')

%% make interpolation objects
[xxGrid,yyGrid,thGrid] = ndgrid(xxIP,yyIP,thetaIP);
fzIP0 = griddedInterpolant(xxGrid,yyGrid,thGrid,permute(zIP,[2,1,3]),'linear','none');
fy0 = griddedInterpolant(xxGrid,yyGrid,thGrid,imag(permute(zIP,[2,1,3])),'linear','none');
fJInv0 = griddedInterpolant(xxGrid,yyGrid,thGrid,abs(permute(f_zIP,[2,1,3])).^(-2) ,'linear','none');
tIP = omega.\asin(thetaIP/thetaMax);
f_tIP2 = cat(3,zeros(size(zzIP)),diff(zIP,1,3)./diff(tIP,1,3),zeros(size(zzIP))); % NB! numerical time differentiation. Find a way to validate!
thetaIP2 = cat(3,thetaIP(1),.5*(thetaIP(2:end)+thetaIP(1:end-1)),thetaIP(end));
f_zIP2 = cat(3,f_zIP(:,:,1),.5*(f_zIP(:,:,2:end)+f_zIP(:,:,1:end-1)),f_zIP(:,:,end));

[xxGrid2,yyGrid2,thGrid2] = ndgrid(xxIP,yyIP,thetaIP2);
ft__fzIP = f_tIP2./f_zIP2;
ft__fzIP(f_tIP2==0) = 0;
ft__fz0 = griddedInterpolant(xxGrid2,yyGrid2,thGrid2,permute(f_tIP2./f_zIP2,[2,1,3]),'linear','none');


% t_ = 0:.01:2*T;
% figure, plot(t_,abs(mod(t_*omega,pi/2)).*sign(cos(t_*omega)),T/4*[1,1],[-1,1]*pi/4,'k')
% figure, plot(t_, interp1(T/4*[-1,1],pi/2*[-1,1],mod(t_+T/4,T/2+1e-12)-T/4).*sign(cos(omega*t_)),  T/4*[1,1].*(1:4)',[-1,1]*pi/4,'k')
% figure, plot(t_, pi*(mod(t_+T/4,T/2+1e-12)-T/4).*sign(cos(omega*t_)),  T/4*[1,1].*(1:4)',[-1,1]*pi/4,'k')
% figure, plot(t_, thetaMax*sin(pi*(mod(t_+T/4,T/2+1e-12)-T/4).*sign(cos(omega*t_))),  T/4*[1,1].*(1:4)',[-1,1]*pi/4,'k')

complex2grid = @(zz,t,ff) ff(real(zz)+0*t,imag(zz)+0*t, thetaMax*sin(omega*t)+0*zz);
fzIP = @(zz,t) complex2grid(zz,t,fzIP0);
map.fy = @(zz,t) complex2grid(zz,t,fy0);
map.fJInv = @(zz,t) complex2grid(zz,t,fJInv0);
map.ft__fz = @(zz,t) complex2grid(zz,t,ft__fz0);


% figure,plot(t_,map.fy(1+.2i,t_))

%% test interpolation with plot
if ~isempty(flapPlotTimes)
    xxEnd = 2*h;
    zz_r = linspace(0,xxEnd,30 )+1i*linspace(-h,wbOverWater,200).';
    zz_b = linspace(0,xxEnd,200)+1i*linspace(-h,wbOverWater,10 ).';

    figure('color','w');hold on; grid on
    plot(0,-wbl,'ok','markersize',8,'linewidth',2);
    axis([-(wbl+wbOverWater)*sin(thetaMax),xxEnd,-h,wbOverWater])
    for iPlot = 1:length(flapPlotTimes), tPlot=flapPlotTimes(iPlot);
        if iPlot>1, delete(hl); end
        hl = [plot(fzIP(zz_r,tPlot),'r','linewidth',1)  
              plot(fzIP(zz_b,tPlot).','b','linewidth',1)
              plot(fzIP(linspace(0,xxEnd,30 ),tPlot).','.-k','linewidth',1)
              plot([0,0,-(wbl+wbOverWater)*tan(thetaMax*sin(omega*tPlot))],[-h,-wbl,wbOverWater],'-k','linewidth',2)];
        drawnow
    end
end

%% tests:
assert(abs(map.fy(1+.2i,.1*T)-map.fy(1+.2i,1.1*T))<1e-12 && abs(map.fy(1+.2i,.1*T)-map.fy(1+.2i,-.9*T))<1e-12)
assert(abs(map.fy(1+.2i,0)-map.fy(1+.2i,T))<1e-12 && abs(map.fy(1+.2i,0)-map.fy(1+.2i,-T))<1e-12)
assert(abs(map.ft__fz(1+.2i,.1*T)-map.ft__fz(1+.2i,1.1*T))<1e-12 && abs(map.ft__fz(1+.2i,.1*T)-map.ft__fz(1+.2i,-.9*T))<1e-12)
assert(abs(map.ft__fz(1+.2i,0)-map.ft__fz(1+.2i,T))<1e-12 && abs(map.ft__fz(1+.2i,0)-map.ft__fz(1+.2i,-T))<1e-12)


% NB! may need altering
map.xi = x;
map.H = h;


% finvIp = scatteredInterpolant(  real(zIP(:)), imag(zIP(:)) , zzIP(:),'linear','none');
% zzS0 = finvIp(x, h0 );
% % interpolate onto regulart xi grid.
% map.xi = linspace(real(zzS0(1)),real(zzS0(nx)),nx)';
% map.H = pi;
% eta0_xiReg = interp1(real(zzS0),imag(zzS0),map.xi);
% xS_xiReg = interp2(real(zzIP),imag(zzIP),real(zIP),map.xi,eta0_xiReg);
% varphiS0 = interp1( [x-L;x;x+L],[phiS0;phiS0;phiS0],xS_xiReg );
    
[varphiS0,eta0_xiReg] = deal(zeros(size(x)));


% return

dim.L  = (map.xi(2)-map.xi(1))*nx/(2*pi);
dim.t = sqrt(dim.L/g);
dim.phi = sqrt(dim.L^3*g);
dim.U = dim.L/dim.t;
param.dim = dim;
param.map = map;

%% Run simulation
tic
if RK4dt~=0
    [t,y] = RK4(@(t,Y) HOS_Taylor(t,Y,param) ,[t0,RK4dt,param.t_end]/dim.t,[varphiS0/dim.phi;eta0_xiReg/dim.L]);
else
    [t,y] = ode45(@(t,Y) HOS_Taylor(t,Y,param) ,[t0,param.t_end]/dim.t,[varphiS0/dim.phi;eta0_xiReg/dim.L],ODEoptions);
end
fprintf('CPU time: %gs\n',toc);
t = t*dim.t;
varphiS = y(:,1:nx)*dim.phi; eta = y(:,nx+1:2*nx)*dim.L;

iNaN = find(isnan(varphiS(:,1)),1,'first');
if ~isempty(iNaN), t(iNaN:end)=[]; varphiS(iNaN:end,:)=[]; eta(iNaN:end,:)=[]; end
clear y



% t_ip = (0:dt:param.t_end);
t_ip = linspace(0,t(end),10);
% t_ip = linspace(0,.9*t(end),10);

nPannel = length(t_ip);
varphiS_ip = interp1(t,varphiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

zS_ip = fzIP(map.xi+1i*eta_ip,t_ip);


[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814]),[.075,.04,.05,.05],[.0,0]); %,'name',sprintf('Conformal; Tramp%g ka=%.3g',TRamp,ka)
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
maxh = max(real(zS_ip(:)));minh = min(real(zS_ip(:)));
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),zS_ip(:,i),'k');
    ylabel(ha(i),sprintf('t = %.2fs',t_ip(i)))
    grid(ha(i),'on');
end
% axis(ha,'equal','tight')
set(ha,'XLim',[minh,maxh],'YLim',[min(imag(zS_ip(:))),max(imag(zS_ip(:)))])
% set(ha,'DataAspectRatio',[1,1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)

% estimate of reflected wave
iRefPan = 9;
ii = real(zS_ip(:,iRefPan))>-.25*L & real(zS_ip(:,iRefPan))<0;
aRef = .5*( max(imag(zS_ip(ii,iRefPan)))-min(imag(zS_ip(ii,iRefPan))));

fileName = sprintf('%s_a1%.4f_aRef%.5f',fileName,ka/k1,aRef); fileName(fileName=='.')='p';
if DO_EXPORT
%     copyfile('./proto_beachRamp.m',[exportPath,'/m/',fileName,'.m']) 
    savefig(hf,[exportPath,'/fig/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],exportFormats{:});
end
if EXPORT_MAT == 1
    wh = whos;
    vars = setdiff({wh.name},{'t','y','varphiS','eta'});
    save([exportPath,'/mat/',fileName],vars{:}); 
elseif EXPORT_MAT == 2
    save([exportPath,'/mat/',fileName]); 
end


%% Volume conserving map
% function [z,df,zz] = fz(xx,yy,theta,h,wbl,wbOverWater)
%     assert(all(diff(xx)>0&&all(diff(yy)>0),'increasing zz vectors assumed');
%     xx = fliplr(xx); % integrating from right-to-left
%     yy = yy+h; 
%     
%     d = wbl+wbOverWater;
%     H = h+wbl+2*wbOverWater;
%     [z,df,zz] = fz0(xx,yy,theta,H,d);
%     
%     z  = fliplr(z ) - 1i*h - real(z(1,end,:));
%     df = fliplr(df);
%     zz = fliplr(zz)-h;
% end
% function [z,df,zz] = fz0(xx,yy,theta,H,d)
% zz = xx + 1i*yy;
% sig0 = H + [-2;-1;1;2]*d;
% 
% sig = zeros(4,1,numel(theta));
% for i = 1:length(theta)
%     sig(:,:,i) = LMFnlsq(@(sig) fixHingeDepth(sig,theta(i),H,d),sig0); % {LMFnlsq,newtonraphson,fminsearch}
% end
% 
% zz2 = zz.^2;
% df = ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);
% df_yh = .5*(df(1:end-1,1,:)+df(2:end,1,:));
% z1 = cumsum([zeros(size(theta));df_yh.*diff(1i*yy)],1);
% df_xh = .5*(df(:,1:end-1,:)+df(:,2:end,:));
% z =  cumsum([z1,df_xh.*diff(xx)],2) ;
% end
% function err = fixHingeDepth(sig,theta,H,d)
% d_xi = 1e-6;
% df_zz2 = @(zz2) ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);
% dy = @(yy) real( df_zz2((1i*yy+d_xi).^2) ); % real beacuse it is integrated by dzz = 1i*dyy;
% y = zeros(5,1); sig0 = [0;sig];
% for i = 1:4
%     y(i+1) = y(i) + integral( dy,sig0(i),sig0(i+1));
% end
% err = y(2:5) - (H+[-2;-1;1;2]*d);
% end



%% flat yy=0 map
% function [z,df,zz] = fz(xx,yy,theta,h,wbl,wbOverWater)
% assert(all(diff(xx)>0&&all(diff(yy)>0),'increasing zz vectors assumed');
% xx = fliplr(xx); % integrating from right-to-left
% 
% thp = theta/pi;
% % d = sqrt(pi).*(wbl+wbOverWater).*sec(theta)./(gamma(1+thp)*gamma(.5-thp));
% zz = xx + 1i*yy;
% % df = (1+d^2./zz.^2).^thp; % k=1
% 
% % stretched variables
% h_   = h   + wbOverWater;
% wbl_ = wbl + wbOverWater;
% zz_  = zz  - 1i*wbOverWater;
% 
% % d = fixHingeDepth(D/H,theta)*H;
% d = 0*theta;
% for i = 1:length(theta)
%     d(i) = fzero(@(d__h) fixHingeDepth_fzero(d__h,wbl_/h_,theta(i)),wbl_/h_)*h_;
% end
% % figure, plot(theta(:),d(:)-linspace(d(1),d(end),length(d))','.-')
% 
% df = (1+csch(pi*zz_./(2*h_)).^2.*sin(pi*d./(2*h_)).^2).^thp; 
% df_yh = .5*(df(1:end-1,1,:)+df(2:end,1,:));
% z1 = cumsum([zeros(size(theta));df_yh.*diff(1i*yy)],1);
% z1 = z1-z1(yy==0);% + 1i*wbOverWater;
% df_xh = .5*(df(:,1:end-1,:)+df(:,2:end,:));
% z =  cumsum([z1,df_xh.*diff(xx)],2) ;
% z = z-real(z(1,end,:));
% 
% z  = fliplr(z );
% df = fliplr(df);
% zz = fliplr(zz);
% 
% end
% 
% function err = fixHingeDepth_fzero(d__h,wbl__h,theta)
% % there's a singularity at zz=-1i*d if theta < 0
% % either stop a yy=-d-delta_singularity
% % or shift integration path d_xi to the right (into the domain):
% % delta_singularity = 1e-6*(theta<0); d_xi = 0;
% delta_singularity = 0; d_xi = 1e-6*(theta<0);
% 
% % nYInt = 10000;
% % yi = linspace(-1,-d__h-delta_singularity,nYInt)';
% % dL = real( (1-csc(pi/2*(yi-1i*d_xi)).^2.*sin(pi/2*d__h).^2).^(theta/pi) );
% % L = sum(.5*(dL(1:end-1)+dL(2:end)) .* diff(yi) )+delta_singularity;
% 
% L = integral(@(y) real((1-csc(pi/2*(y-1i*d_xi)).^2.*sin(pi/2*d__h).^2).^(theta/pi)),-1,-d__h-delta_singularity)+delta_singularity;
% err = wbl__h-(1-L);
% end
















% % function d__h = fixHingeDepth(wbl__h,theta)
% % maxIt = 100;
% % nYInt = 500;
% % d__h = wbl__h+0*theta;
% % 
% % % there's a singularity at zz=-1i*d if theta < 0
% % % either stop a yy=-d-delta_singularity
% % % or shift integration path d_xi to the right (into the domain):
% % % delta_singularity = 1e-6*(theta<0); d_xi = 0;
% % delta_singularity = 0; d_xi = 1e-6*(theta<0);
% % 
% % for i = 1:maxIt
% %     
% %     dy = (-d__h-delta_singularity+1)/(nYInt-1);
% %     yi = -1 + cumsum((0:nYInt-1)'.*dy,2);
% % 
% %     dL = real( (1-csc(pi/2*(yi-1i*d_xi)).^2.*sin(pi/2*d__h).^2).^(theta/pi) ); 
% %     L = sum(.5*(dL(1:end-1,:,:)+dL(2:end,:,:)) .* diff(yi,1,1), 1 )+ delta_singularity ;
% %     fprintf('i %d: L = %g, 1-wbl/h=%g\n',i,max(L(:)),1-wbl__h)
% % 
% %     delta_d__h = wbl__h-(1-L);
% %     
% %     iChange = abs(delta_d__h)>1e-9;
% %     d__h(iChange) = d__h(iChange) + delta_d__h(iChange);
% %     if ~any(iChange), return; end
% % end
% % 
% % [~,iMax] = max(abs(delta_d__h(:)));
% % warning('failed to find solution withing %d iterations. Biggest last error %g at theta = %.2g deg.',maxIt,delta_d__h(iMax),theta(iMax)*180/pi)
% % end
% 
% 
% % 
% % function fInv_t = getFInv_t(xx,yy,theta_t,h,wbl,wbOverWater)
% % 
% % xx = fliplr(xx); % integrating from right-to-left
% % yy = flipud(yy); % integrating from right-to-left
% % zz = xx + 1i*yy;
% % H = h + wbOverWater;
% % dRep = wbl + wbOverWater;
% % dg = -log(1+csch(pi*zz./(2*H)).^2.*sin(pi*dRep./(2*H)).^2);% time derivative term
% % dg_yh = .5*(dg(1:end-1,1)+dg(2:end,1));
% % g1 = cumsum([0;dg_yh.*diff(1i*yy)],1);
% % dg_xh = .5*(dg(:,1:end-1)+dg(:,2:end));
% % g =  cumsum([g1,dg_xh.*diff(xx)],2) ;
% % g = g-g(end,end);
% % fInv_t = rot90(g,2).*theta_t/pi; 
% % end
