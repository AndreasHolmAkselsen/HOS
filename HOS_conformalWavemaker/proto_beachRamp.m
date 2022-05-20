clear 
global timeReached 
timeReached = 0; 

g = 9.81;
% g = 1.0;

%% input

%%%%%%%%%%%%%%%%%%%% 5m depth, 1m deep extension
% % 45 deg 
% xx_b = [0,nan];    % xi-coordinate of "begining of" edge
% H = [5,2.4,1]; % plateau levels
% theta = [1,.5]*pi/2; % slope angles (positive values)

% % 30 deg
% xx_b = [0,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,2.25,1]; % plateau levels
% theta = [90,30]*pi/180; % slope angles (positive values)

% % straight step
% xx_b = [0];    % xi-coordinate of "begining of" edge
% H = [5,1]; % plateau levels
% theta = [1]*pi/2; % slope angles (positive values)
%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%% 10m depth, 2m deep extension sticking out
% % 30 deg
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [10,9.1,4.6,2]; % plateau levels
% theta = [90,180,30]*pi/180; % slope angles (positive values)

% % 45 deg
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [10,8.3,4.9,2]; % plateau levels
% theta = [90,180,45]*pi/180; % slope angles (positive values)

% % straight step
% xx_b = [0];    % xi-coordinate of "begining of" edge
% H = [10,2]; % plateau levels
% theta = 90*pi/180; % slope angles (positive values)
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%  beach tests, 5m depth, 2m deep extension sticking out
% % 30 deg
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,4.95,4.33,2]; % plateau levels
% theta = [90,180,30]*pi/180; % slope angles (positive values)

% % 45 deg
% xx_b = [0,nan,nan];  %2.7726  % xi-coordinate of "begining of" edge
% H = [5,4.95,4.43,2]; % plateau levels
% theta = [90,180,45]*pi/180; % slope angles (positive values)

% 90 deg
xx_b = [0];  %2.7726  % xi-coordinate of "begining of" edge
H = [5,2]; % plateau levels
theta = [90]*pi/180; % slope angles (positive values)
%%%%%%%%%%%%


% interpolation grid for conformal map
nArrX_far = 1500;   % # horizontal points (far field)
nArrX_near = 500;  % # horizontal points (near field)
nArrYDown = 100;    % # vertical points

H_IC = H(1);


% for wave breaking example
param.DO_PADDING = 0;
RK4dt = 0;%2.5e-3; % set to zero to use ODE45
ka = .05; % linear wave steepness
%  and  NT_dt =  5 /9/T; lambda = 2*pi; g=1;


relTolODE = 1e-4;% 1e-8;
N_SSGW = 2^12; % number of modes in SSGW solution

% Plot & export options
DO_EXPORT = 1;
EXPORT_MAT = 1;
PLOT_MAP = 0;
PLOT_INTERPOLATION_MAP = 0;
exportPrefix = 'beach_';
exportPath = './figures/';
exportFormatsMap = {'-pdf','-png'};
exportFormats = {'-png','-pdf','-m2'};


% Wave init specification
% specifying right-side relative depth
% kH_R = .66;
% k_R = kH_R/H(end);
% omega = sqrt(g*k_R*tanh(kH_R)); T = 2*pi/omega;
% k1 = findWaveNumbers(omega,H(1),0,0);

% % specifying period
T = 3.5; omega=2*pi/T; k1 = findWaveNumbers(omega,H(1),0,0);


NWaves = 60;
lambda1 = 2*pi/k1;
L = lambda1*NWaves;
xLR = [-L/2,L/2];
% width_x = width_x__L*L;


% T = 1; omega=2*pi/T;
%  if U_curr==0, k0=omega^2/g;else, k0=(g+2*U_curr*omega-sqrt(g^2+4*g*U_curr*omega))/(2*U_curr^2);end;lambda=2*pi/k0;


% Simulation/plotting time
NT_dt = 5;
dt = NT_dt*T;
param.t_end = 9*dt;
    
% Initial conditions
INIT_WAVE_TYPE = 'SSGW';  % 'SSGW' or 'linear'
packageWidth__L = .05;  % set to inf if not simulating wave packets
packageCentre__L = -.2;
% packageCentre__L = -1/3;

nx__wave = 2^6;
param.M = 5; 
TRamp = 0*T;
param.nonLinRamp = @(t) max(0,1-exp(-(t/TRamp)^2));
% param.iModeCut = 2*(NWaves + (param.M+5));

param.iModeCut = inf;
param.kd__kmax = .5;
param.rDamping = .25;



%% Simulation
nx = nx__wave*NWaves;

dx = L/nx;
% x = (0:nx-1)'*dx;
x = linspace(xLR(1),xLR(2),nx+1)';x(end)=[];
fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-param.iModeCut/ (nx/2),0) )
packet = exp(-((x/L-packageCentre__L)/packageWidth__L).^2);

t0 = 0;
initialStepODE = 1e-3*T;
xk1 = k1.*x;
phaseAng = 0*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);%,'MaxStep',1e-4);

Hstr = sprintf('%.2f_',H); thetaStr = sprintf('%.0f_',theta*180/pi);
fileName = sprintf('%s%s_T%.2f_ka%.2g_M%d_H%stheta%sNw%d_dt%.3gT_nx%d_pad%d_ikCut%.4g_Md%.2g_r%.2g',exportPrefix,INIT_WAVE_TYPE,T,ka,param.M,Hstr,thetaStr,NWaves,NT_dt,nx,param.DO_PADDING,param.iModeCut,param.kd__kmax,param.rDamping); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./proto_beachRamp.m',[exportPath,'/m/',fileName,'.m']) 
end


%% Initial conditions:
switch INIT_WAVE_TYPE
    case 'linear'
        h0 = ka/k1*(cos(xk1-phaseAng));
        phiS0 = ka/k1.*g/omega*sin(xk1-phaseAng);
        
    case 'SSGW'
        [zIC,dwdz,PP] = SSGW(k1*H_IC,ka,N_SSGW);
        
        if isinf(PP(1)), L_scale = 1/k0; else, L_scale = H_IC; end
        out.c_e = PP(4)*sqrt(g*L_scale); % phase velocity observed from where the meam velocity at the bed is zero
        out.c_s = PP(5)*sqrt(g*L_scale); % mean flow velocity (phase velocity in frame without mean flow)
        out.k = PP(2)/L_scale;
        zIC = zIC*L_scale;
        
%         % to move wave to centre (optional)
%         z = [ z(N_SSGW+1:end)-lambda/2 ; z(1:N_SSGW)+lambda/2 ];
%         dwdz = [ dwdz(N_SSGW+1:end); dwdz(1:N_SSGW) ];
        
        % duplicate across domain.
        zIC = reshape(repmat(zIC,1,NWaves)+lambda1*(0:NWaves-1),[],1) + x(1);
        dwdz = repmat(dwdz,NWaves,1);
        dwdz = dwdz*sqrt(g*L_scale);
        
        n = 2*N_SSGW*NWaves;
        z_m = .5*(zIC(1:n-1)+zIC(2:n));
        dwdz0_m = .5*(dwdz(1:n-1)+dwdz(2:n))+out.c_e;
        wIC = [0;cumsum( dwdz0_m.*diff(zIC))];
        wIC = wIC-mean(wIC);
        
        % if z(1)<2*eps&&z(1)>-2*eps, z(1)=1i*imag(z(1));end
        zIC = [zIC(end)-L;zIC;zIC(1)+L]; % extend with ghost nodes
        wIC = [wIC(end);wIC;wIC(end)];
        h0 = interp1(real(zIC),imag(zIC),x,'linear',nan);
        phiS0 = interp1(real(zIC),real(wIC),x,'linear',nan);
        fft_h = fftshift(fft(h0));
        if sum(abs(fft_h(1:floor(end/4))))>.01*sum(abs(fft_h))
            warning('Initial conition may not have been found. Verify that solution exists.')
        end
%         phiS00 = ka/k0.*g/omega*sin(xk0-phaseAng);
%         h00 = ka/k0*(cos(xk0-phaseAng));
%         % phi = ka/k0.*g/omega*sin(xk0-phaseAng)*cosh(k*(h+z))/cosh(k*h);
%         u0 = ka/k0.*g/omega*k0*cos(xk0-phaseAng);
%         v0 =  ka/k0.*g/omega*k0*sin(xk0-phaseAng)*tanh(k0*H_IC);
%         figure('color','w')
%         subplot(311), plot(x,h0,'-',x,h00,'--');ylabel('\eta'); grid on; legend('SSGW','linear')
%         subplot(312), plot(x,phiS0,'-',x,phiS00,'--');ylabel('\phi^S'); grid on; legend('SSGW','linear')
%         subplot(313), plot(x,u0,'-r',x,v0,'-b',real(zIC),real(dwdz)+out.c_e,'--r',real(zIC),-imag(dwdz),'--b');ylabel('velocity'); grid on
        
    otherwise % input data file assumed
        filePath = ['./IC/',INIT_WAVE_TYPE,'.mat'];
        assert(isfile(filePath));
        load(filePath)
        assert(x0(end)<=x(end))
        h0 = [h0;zeros(nx-length(x0),1)];
        varphiS0 = [varphiS0;zeros(nx-length(x0),1)];
%         figure, plot(x,eta,x,phiS,'--')
end
% Adjust in case of wave packets.
h0 = h0.*packet;
phiS0 = phiS0.*packet;

% dim.L  = L/(2*pi);
% dim.t = sqrt(dim.L/g);
% dim.phi = sqrt(dim.L^3*g);
% dim.U = dim.L/dim.t;


%%  map preparation

% Create inverse function through scattered interpolation:

for iNan = find(isnan(xx_b))
    xx_b(iNan) = xx_b(iNan-1) - pi./theta(iNan).*log(H(iNan+1)./H(iNan));
end
zzRoots = [xx_b-1i*pi;xx_b+pi./theta.*log(H(2:end)./H(1:end-1))-1i*pi].';
crudeScale = 1.5*pi/min(H);
yyUpper = 2*max(h0)*crudeScale;
% xxIP = linspace(xLR(1),xLR(2),nArrX)*crudeScale; 

xxIP_far = linspace(xLR(1),xLR(2),nArrX_far)*crudeScale; 
% xxIP_near = linspace(xx_b(1)-delta_xx_nearField,xx_b(end)+delta_xx_nearField,nArrX_near);
% cSq = max(H(2:end)./H(1:end-1),H(1:end-1)./H(2:end)).^(pi/theta);
% xxIP_near = linspace(xx_b(1)-1.5*log(cSq(1)-1),xx_b(end)+1.5*log(cSq(end)-1),nArrX_near);

% singularities are in the zz-plane located at
% xi_j - 1i*pi and xi_j + log(c^2) - 1i*pi

% log_cSq = pi./theta.*log(H(2:end)./H(1:end-1));
% xxIP_near = linspace( xx_b(1)+min(0,log_cSq(1))-.5*abs(log_cSq(1)), xx_b(end)+max(0,log_cSq(end))+.5*abs(log_cSq(end)), nArrX_near);

xxRoots = real(zzRoots);
% xxIP_near = linspace( min(real(xxRoots(1,:)))-20*abs(diff(xxRoots(1,:))), max(real(xxRoots(end,:)))+5*abs(diff(xxRoots(end,:))), nArrX_near);
dxxRoot = max(xxRoots(:))-min(xxRoots(:));
xxIP_near = linspace( min(real(xxRoots(1,:)))-5*dxxRoot, max(real(xxRoots(end,:)))+5*dxxRoot, nArrX_near);


xxIP = [xxIP_far(xxIP_far<xxIP_near(1)), xxIP_near, xxIP_far(xxIP_far>xxIP_near(end))];
H = shiftdim(H,-1); theta = shiftdim(theta,-1);xx_b = shiftdim(xx_b,-1);

% yyIP = linspace(-pi,1.4*max(h0)*crudeScale,nArrYDown)';
yyIP = linspace(-pi,0,nArrYDown)'; % ensure that we capture the line yy=0
dy = yyIP(2)-yyIP(1);
yyIP = [yyIP; (dy:dy:yyUpper)' ];
assert(yyIP(nArrYDown)==0)

[zzIP,dfIP,zIP] = fz(xxIP,yyIP,H,theta,xx_b);
assert(real(zIP(nArrYDown,1))<x(1)&&real(zIP(nArrYDown,end))>x(end),'Physical domain [%.3g,%.3g] out of range of interpolation range [%.3g,%.3g]. Extend interpolation range.',x(1),x(end),real(zIP(nArrYDown,1)),real(zIP(nArrYDown,end)))

% plot interpolation basis to inspect resolution:
if PLOT_INTERPOLATION_MAP
    figure('color','w');hold on
    plot(zIP(:,1:round(end/30):end),'r','linewidth',1);  plot(zIP(1:round(end/10):end,:).' ,'b')
    minIz = min(imag(zIP(1,:))); patch(real(zIP(1,[1,1:end,end])),[1.1*minIz,imag(zIP(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
end


if ~all(abs(imag(zIP(nArrYDown,:)))<1e-3*max(H))
    warning('Map appears not to have a flat surface at yy=0. max|y(yy=0)| = %g',max(abs(imag(zIP(nArrYDown,:)))))
end
% assert(all(abs(imag(zIP(nArrYDown,:)))<1e-3*max(H)),'Map appears not to have a flat surface at yy=0')
xxLR = interp1(real(zIP(nArrYDown,:)),xxIP,xLR);

% iTrim = xxIP>=xxLR(1)&xxIP<=xxLR(2);
iTrim = (find(xxIP<xxLR(1),1,'last')-1):(find(xxIP>xxLR(2),1,'first')+1);
zzIP = zzIP(:,iTrim); zIP = zIP(:,iTrim); dfIP = dfIP(:,iTrim); xxIP = xxIP(:,iTrim);
assert(all(real(zzIP(:,1))<xxLR(1)) && all(real(zzIP(:,end))>xxLR(2)))


fzIP0 = griddedInterpolant(real(zzIP).',imag(zzIP)',zIP.','linear','none');
fy0 = griddedInterpolant(real(zzIP).',imag(zzIP).',imag(zIP).','linear','none');
fJInv0 = griddedInterpolant(real(zzIP).',imag(zzIP).', abs(dfIP).'.^(-2) ,'linear','none');

complex2grid = @(zz,ff) ff( (real(zz)+0*zz).', (imag(zz)+0*zz).').';
fzIP = @(zz) complex2grid(zz,fzIP0);
map.fy = @(zz) complex2grid(zz,fy0);
map.fJInv = @(zz) complex2grid(zz,fJInv0);


finvIp = scatteredInterpolant(  real(zIP(:)), imag(zIP(:)) , zzIP(:),'linear','none');
zzS0 = finvIp(x, h0 );
% interpolate onto regulart xi grid.
map.xi = linspace(real(zzS0(1)),real(zzS0(nx)),nx)';
map.H = pi;
eta0_xiReg = interp1(real(zzS0),imag(zzS0),map.xi);
xS_xiReg = interp2(real(zzIP),imag(zzIP),real(zIP),map.xi,eta0_xiReg);
varphiS0 = interp1( [x-L;x;x+L],[phiS0;phiS0;phiS0],xS_xiReg );
    
% figure, contourf(real(zIP),imag(zIP),map.fJInv(zIP));colorbar
% figure, subplot(121); plot(xS_xiReg,map.fy(map.xi+1i*eta0_xiReg),x,h0,'--'); 
% subplot(122);plot(xS_xiReg,varphiS0,x,phiS0,'--');




if PLOT_MAP
    hf_map = figure('color','w','position',[436 63 637 600]); 

    % plot the z-plane
    haz = subplot(211); hold on;
%     haz = axes; hold on;
    title('z-plane');xlabel('x');ylabel('i y');box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    
%     xxPlot = linspace(xxLR(1),xxLR(2),200);
%     xxPlot = xxIP;
%     xxPlot = xxIP(xxIP>=xxLR(1)&xxIP<=xxLR(2));
    nxxPlot = 1000;
    xxPlot = xxIP(1:round(length(xxIP)/nxxPlot):end);
%     zPhi = fzIP0(repmat(linspace(xxLR(1),xxLR(2),20)',1,200),repmat(linspace(-pi,yyUpperPlot,200),20,1)).';
%     zPsi = fzIP0(repmat(linspace(xxLR(1),xxLR(2),200)',1,10),repmat(linspace(-pi,yyUpperPlot,10),200,1)).';
    zPhi = fzIP(linspace(xxLR(1),xxLR(2),20)+1i*linspace(-pi,yyUpper,200).');
    zPsi = fzIP(xxPlot+1i*linspace(-pi,yyUpper,10).');
    plot(zPhi,'r','linewidth',1); hold on; plot(zPsi.' ,'b')
    minIz = min(imag(zPsi(1,:)));
    patch(real(zPsi(1,[1,1:end,end])),[1.1*minIz,imag(zPsi(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');%,'FaceAlpha',.5
    xlabel('x');ylabel('i y');
    plot(x,h0,'k','linewidth',1.5);
    
    % plot the zz-plane
    hazz = subplot(212); hold on    
    title('\zeta-plane'); xlabel('\xi');ylabel('i \sigma'); box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    zzPlot = linspace(xxLR(1),xxLR(2),100) + 1i*linspace(-pi,yyUpper,100)';
    zPlot = fzIP(zzPlot); 
    contour(real(zzPlot),imag(zzPlot),real(zPlot),20,'r','linewidth',1);
    contour(real(zzPlot),imag(zzPlot),imag(zPlot),10,'b','linewidth',1);
    plot(zzS0,'k','linewidth',1.5);
    
    
    % repeat for a 1-to-1 plot
%     xxLR = [-5,10];
    xxLRnear = xxIP_near([1,end]);
    hf_mapZoom = figure('color','w','position',[436 63 637 600]); 
    % plot the z-plane
    haz = subplot(211); hold on;
%     haz = axes; hold on;
    title('z-plane');xlabel('x');ylabel('i y');box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    zPhi = fzIP(linspace(xxLRnear(1),xxLRnear(2),20)+1i*linspace(-pi,yyUpper,nxxPlot).');
    zPsi = fzIP(linspace(xxLRnear(1),xxLRnear(2),nxxPlot)+1i*linspace(-pi,yyUpper,10).');
    plot(zPhi,'r','linewidth',1); hold on; plot(zPsi.' ,'b','linewidth',1)
    minIz = min(imag(zPsi(1,:)));
    patch(real(zPsi(1,[1,1:end,end])),[1.1*minIz,imag(zPsi(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');%,'FaceAlpha',.5
    xlabel('x');ylabel('i y');
    plot(x,h0,'k','linewidth',1.5);
    plot(fzIP(zzRoots+.025i),'ro'); % mark singularities
    axis(haz,'equal','tight');xlim(haz,real(fzIP(xxLRnear)));

        
    % plot the zz-plane
    hazz = subplot(212); hold on    
    title('\zeta-plane'); xlabel('\xi');ylabel('i \sigma'); box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    zzPlot = linspace(xxLRnear(1),xxLRnear(2),100) + 1i*linspace(-pi,yyUpper,100)';
    zPlot = fzIP(zzPlot); 
    contour(real(zzPlot),imag(zzPlot),real(zPlot),20,'r','linewidth',1);
    contour(real(zzPlot),imag(zzPlot),imag(zPlot),10,'b','linewidth',1);
    plot(zzS0,'k','linewidth',1.5);
    plot(zzRoots,'ro'); % mark singularities
    axis(hazz,'equal');xlim(hazz,xxLRnear);
    

    if DO_EXPORT
        Hstr = sprintf('%.2f_',H); thetaStr = sprintf('%.0f_',theta*180/pi);
        fileNameMap = sprintf('%s%s_ka%.2g_H%stheta%sNw%d',exportPrefix,INIT_WAVE_TYPE,ka,Hstr,thetaStr,NWaves); fileNameMap(fileNameMap=='.')='p';
        export_fig(hf_map,['./figures/map/map_',fileNameMap],exportFormatsMap{:})
        savefig(hf_map,['./figures/fig/map_',fileNameMap])
        export_fig(hf_mapZoom,['./figures/map/mapZoom_',fileNameMap],exportFormatsMap{:})
        savefig(hf_mapZoom,['./figures/fig/mapZoom_',fileNameMap])
    end
    
end


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



% t_ip = (0:dt:param.t_end)';
t_ip = linspace(0,t(end),10).';
% t_ip = linspace(0,.9*t(end),10).';

nPannel = length(t_ip);
varphiS_ip = interp1(t,varphiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

zS_ip = fzIP(map.xi+1i*eta_ip);


[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('Conformal; Tramp%g ka=%.3g',TRamp,ka)),[.075,.04,.05,.05],[.0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
maxh = max(real(zS_ip(:)));minh = min(real(zS_ip(:)));
% zSingular = fzIP([xx_b-1i*pi,xx_b+log_cSq-1i*pi]+.025i);
zSingular = fzIP(zzRoots);
% x_b = real(fzIP(xx_b-1i*pi));
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),zS_ip(:,i),'k');
    ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),param.nonLinRamp(t_ip(i))))
    grid(ha(i),'on');
    plot(ha(i),[1;1].*real(zSingular(:)).',[minh;maxh],'--k'); 
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


function [zz,df,z] = fz(xx,yy,H,theta,xx_b)
    assert(isvector(xx)&&isvector(yy));
    xx = xx(:)'; yy = yy(:);
    assert(yy(2)>yy(1))
    yy = flipud(yy);
    zz = xx + 1i*yy;
    
    % compute df for each corner
    lambda = exp(zz-xx_b);
    df_i = ((lambda+1)./(lambda+(H(2:end)./H(1:end-1)).^(pi/theta))).^(theta/pi); % df ~ 1/tau
    
    xi_cut = 50; % xi_value at which to follow asymptote in step.
    iPlus = real(zz-xx_b) > xi_cut;
    iMinus = real(zz-xx_b) < -xi_cut;
    df_i(iPlus) = 1;
    temp = H(1:end-1)./H(2:end)+0*zz;
    df_i(iMinus) = temp(iMinus);
   
    df = prod(df_i,3);
    [ny,nx] = size(zz);
    df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
    df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
    z = cumsum([0,df_xh.*diff(xx)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*yy)],1);
%     z0 = interp2(xx,yy,real(z),0,0)+1i*interp2(xx,yy,imag(z),0,0);
    z0 = interp2(xx,yy,real(z),0,-pi)+1i*interp2(xx,yy,imag(z),0,0);
    z = z-z0; % orientation constant
    
    if isscalar(xx_b)
        K = H(2)./pi;
    else
        K = 1/(-min(imag(z(ny,:))))*max(H);% scaling constant
    end
    z = K*z;
    df = K*df;
        
    [zz,df,z] = deal(flipud(zz),flipud(df),flipud(z));
end

% Obs! df = prod(df_i) -> ddf = D[df] ~= prod(ddf_i)!!!
% function [zz,ddf] = ddfz(xx,yy)
%     global H theta xx_b
%     assert(isvector(xx)&&isvector(yy));
%     xx = xx(:)'; yy = yy(:);
%     assert(yy(2)>yy(1))
%     yy = flipud(yy);
%     zz = xx + 1i*yy;
%     
%     % compute df for each corner
%     c = (H(2:end)./H(1:end-1)).^(pi/theta/2);
%     lambda = exp(zz-xx_b);
%     tau_i = ((lambda+c.^2)./(lambda+1)).^(theta/pi); % df ~ 1/tau
%     
%     ddf_i = tau_i.*lambda.*(c.^2-1)./(c.^2+lambda).^2;
%     xi_cut = 50; % xi_value at which to follow asymptote in step.
%     ddf_i( real(zz-xx_b) > xi_cut |  real(zz-xx_b) < -xi_cut ) = 0;
%     ddf = prod(ddf_i,3);
%     [zz,ddf] = deal(flipud(zz),flipud(ddf));    
% end
