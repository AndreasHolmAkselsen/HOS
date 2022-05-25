clear 
global timeReached 
timeReached = 0; 

g = 9.81;
% g = 1.0;

%% input



% xx_b = [-50,50];    % xi-coordinate of "begining of" edge
% H = [1,.5,1]; % plateau levels
% theta = [.5,.5]*pi/2; % slope angles (positive values)


xx_b = 0;    % xi-coordinate of "begining of" edge
H = [1,.5]; % plateau levels
% theta = [.1]*pi/2; % slope angles (positive values)
theta = [1]*pi/2; % slope angles (positive values)

% interpolation grid for conformal map
nArrX_far = 1000;   % # horizontal points (far field)
nArrX_near = 2000;  % # horizontal points (near field)
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
PLOT_MAP = 1;
exportPrefix = '';
exportPath = './figures/';
i_detailedPlot = []; %plot contour plots of frame i. Leave empty to skip


% Wave init specification

kH1 = 1;
NWaves = 60;
k0 = kH1/H_IC;
lambda = 2*pi/k0;
L = lambda*NWaves;
xLR = [-L/2,L/2];
% width_x = width_x__L*L;

% omega = (1+.5*ka^2)*sqrt(g*k0*tanh(k0*H));
% omega = k0*U_curr+sqrt(g*k0*tanh(k0*H)); T = 2*pi/omega;
% omega = sqrt(g*k0*tanh(k0*.5*(H1+H2))); T = 2*pi/omega;
omega = sqrt(g*k0*tanh(k0*(H_IC))); T = 2*pi/omega;
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
xk0 = k0.*x;
phaseAng = 0*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);%,'MaxStep',1e-4);



%% Initial conditions:
switch INIT_WAVE_TYPE
    case 'linear'
        h0 = ka/k0*(cos(xk0-phaseAng));
        phiS0 = ka/k0.*g/omega*sin(xk0-phaseAng);
        
    case 'SSGW'
        [zIC,dwdz,PP] = SSGW(k0*H_IC,ka,N_SSGW);
        
        if isinf(PP(1)), L_scale = 1/k0; else, L_scale = H_IC; end
        out.c_e = PP(4)*sqrt(g*L_scale); % phase velocity observed from where the meam velocity at the bed is zero
        out.c_s = PP(5)*sqrt(g*L_scale); % mean flow velocity (phase velocity in frame without mean flow)
        out.k = PP(2)/L_scale;
        zIC = zIC*L_scale;
        
%         % to move wave to centre (optional)
%         z = [ z(N_SSGW+1:end)-lambda/2 ; z(1:N_SSGW)+lambda/2 ];
%         dwdz = [ dwdz(N_SSGW+1:end); dwdz(1:N_SSGW) ];
        
        % duplicate across domain.
        zIC = reshape(repmat(zIC,1,NWaves)+lambda*(0:NWaves-1),[],1) + x(1);
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


crudeScale = 1.5*pi/min(H);
yyUpper = 2*max(h0)*crudeScale;
% xxIP = linspace(xLR(1),xLR(2),nArrX)*crudeScale; 

xxIP_far = linspace(xLR(1),xLR(2),nArrX_far)*crudeScale; 
% xxIP_near = linspace(xx_b(1)-delta_xx_nearField,xx_b(end)+delta_xx_nearField,nArrX_near);
% cSq = max(H(2:end)./H(1:end-1),H(1:end-1)./H(2:end)).^(pi/theta);
% xxIP_near = linspace(xx_b(1)-1.5*log(cSq(1)-1),xx_b(end)+1.5*log(cSq(end)-1),nArrX_near);

% singularities are in the zz-plane located at
% xi_j - 1i*pi and xi_j + log(c^2) - 1i*pi

log_cSq = pi./theta.*log(H(2:end)./H(1:end-1));
xxIP_near = linspace( xx_b(1)+min(0,log_cSq(1))-.5*abs(log_cSq(1)), xx_b(end)+max(0,log_cSq(end))+.5*abs(log_cSq(end)), nArrX_near);

xxIP = [xxIP_far(xxIP_far<xxIP_near(1)), xxIP_near, xxIP_far(xxIP_far>xxIP_near(end))];
H = shiftdim(H,-1); theta = shiftdim(theta,-1);xx_b = shiftdim(xx_b,-1);

% yyIP = linspace(-pi,1.4*max(h0)*crudeScale,nArrYDown)';
yyIP = linspace(-pi,0,nArrYDown)'; % ensure that we capture the line yy=0
dy = yyIP(2)-yyIP(1);
yyIP = [yyIP; (dy:dy:yyUpper)' ];
assert(yyIP(nArrYDown)==0)

[zzIP,dfIP,zIP] = fz(xxIP,yyIP,H,theta,xx_b);
assert(real(zIP(nArrYDown,1))<x(1)&&real(zIP(nArrYDown,end))>x(end),'Physical domain [%.3g,%.3g] out of range of interpolation range [%.3g,%.3g]. Extend interpolation range.',x(1),x(end),real(zIP(nArrYDown,1)),real(zIP(nArrYDown,end)))

% hfz = figure('color','w');hold on
% contour(real(zIP),imag(zIP),real(zzIP),10,'r','linewidth',1);
% contour(real(zIP),imag(zIP),imag(zzIP),10,'b','linewidth',1);
% minIz = min(imag(zIP(1,:)));
% hp = patch(real(zIP(1,[1,1:end,end])),[1.1*minIz,imag(zIP(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');

if ~all(abs(imag(zIP(nArrYDown,:)))<1e-3*max(H))
    warning('Map appears not to have a flat surface at yy=0')
end
% assert(all(abs(imag(zIP(nArrYDown,:)))<1e-3*max(H)),'Map appears not to have a flat surface at yy=0')
xxLR = interp1(real(zIP(nArrYDown,:)),xxIP,xLR);

fzIP0 = griddedInterpolant(real(zzIP).',imag(zzIP)',zIP.','linear','none');
fy0 = griddedInterpolant(real(zzIP).',imag(zzIP).',imag(zIP).','linear','none');
fJInv0 = griddedInterpolant(real(zzIP).',imag(zzIP).', abs(dfIP).'.^(-2) ,'linear','none');

complex2grid = @(zz,ff) ff( (real(zz)+0*zz).', (imag(zz)+0*zz).').';
fzIP = @(zz) complex2grid(zz,fzIP0);
map.fy = @(zz) complex2grid(zz,fy0);
map.fJInv = @(zz) complex2grid(zz,fJInv0);
assert(all(zzIP(:,1)<xLR(1)) && all(zzIP(:,end)>xLR(2)))


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
    title('z-plane');xlabel('x');ylabel('i y');box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    
%     zPhi = fzIP0(repmat(linspace(xxLR(1),xxLR(2),20)',1,200),repmat(linspace(-pi,yyUpperPlot,200),20,1)).';
%     zPsi = fzIP0(repmat(linspace(xxLR(1),xxLR(2),200)',1,10),repmat(linspace(-pi,yyUpperPlot,10),200,1)).';
    zPhi = fzIP(linspace(xxLR(1),xxLR(2),20)+1i*linspace(-pi,yyUpper,200).');
    zPsi = fzIP(linspace(xxLR(1),xxLR(2),200)+1i*linspace(-pi,yyUpper,10).');
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
    xxLR = xxIP_near([1,end]);
    hf_mapZoom = figure('color','w','position',[436 63 637 600]); 
    % plot the z-plane
    haz = subplot(211); hold on;
    title('z-plane');xlabel('x');ylabel('i y');box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    zPhi = fzIP(linspace(xxLR(1),xxLR(2),20)+1i*linspace(-pi,yyUpper,200).');
    zPsi = fzIP(linspace(xxLR(1),xxLR(2),200)+1i*linspace(-pi,yyUpper,10).');
    plot(zPhi,'r','linewidth',1); hold on; plot(zPsi.' ,'b')
    minIz = min(imag(zPsi(1,:)));
    patch(real(zPsi(1,[1,1:end,end])),[1.1*minIz,imag(zPsi(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');%,'FaceAlpha',.5
    xlabel('x');ylabel('i y');
    plot(x,h0,'k','linewidth',1.5);
    zzSingular = [xx_b-1i*pi,xx_b+log_cSq-1i*pi];
    plot(fzIP(zzSingular+.025i),'ro'); % mark singularities
    % plot the zz-plane
    hazz = subplot(212); hold on    
    title('\zeta-plane'); xlabel('\xi');ylabel('i \sigma'); box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    zzPlot = linspace(xxLR(1),xxLR(2),100) + 1i*linspace(-pi,yyUpper,100)';
    zPlot = fzIP(zzPlot); 
    contour(real(zzPlot),imag(zzPlot),real(zPlot),20,'r','linewidth',1);
    contour(real(zzPlot),imag(zzPlot),imag(zPlot),10,'b','linewidth',1);
    plot(zzS0,'k','linewidth',1.5);
    plot(zzSingular,'ro'); % mark singularities
    
    axis(haz,'equal','tight');xlim(haz,real([zPlot(1,1),zPlot(1,end)]));
    axis(hazz,'equal');xlim(hazz,xxLR);
    

    if DO_EXPORT
        fileNameMap = sprintf('%s%s_ka%.2g_H%.2f_%.2f_nH%d_ang1_%.2g_Nw%d',exportPrefix,INIT_WAVE_TYPE,ka,H(1),H(2),length(H),2*theta/pi,NWaves); fileNameMap(fileNameMap=='.')='p';
        export_fig(hf_map,['./figures/map/map_',fileNameMap],'-png','-pdf','-m2')
        savefig(hf_map,['./figures/fig/map_',fileNameMap])
        export_fig(hf_mapZoom,['./figures/map/mapZoom_',fileNameMap],'-png','-pdf','-m2')
        savefig(hf_mapZoom,['./figures/fig/mapZoom_',fileNameMap])
    end
    
end

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
zSingular = fzIP([xx_b-1i*pi,xx_b+log_cSq-1i*pi]+.025i);
% x_b = real(fzIP(xx_b-1i*pi));
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),zS_ip(:,i),'k');
    ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),param.nonLinRamp(t_ip(i))))
    grid(ha(i),'on');
    plot(ha(i),[1;1].*real(zSingular),[minh;maxh],'--k'); 
end
% axis(ha,'equal','tight')
set(ha,'XLim',[minh,maxh],'YLim',[min(imag(zS_ip(:))),max(imag(zS_ip(:)))])
% set(ha,'DataAspectRatio',[1,1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)

fileName = sprintf('%s%s_ka%.2g_M%d_H%.2f_%.2f_nH%d_ang1_%.2g_Nw%d_dt%.3gT_nx%d_pad%d_ikCut%.4g_Md%.2g_r%.2g',exportPrefix,INIT_WAVE_TYPE,ka,param.M,H(1),H(2),length(H),2*theta/pi,NWaves,NT_dt,nx,param.DO_PADDING,param.iModeCut,param.kd__kmax,param.rDamping); fileName(fileName=='.')='p';

if DO_EXPORT
    copyfile('./proto_sloping.m',[exportPath,'/m/',fileName,'.m']) 
    savefig(hf,[exportPath,'/fig/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf','-png');
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
    c = (H(2:end)./H(1:end-1)).^(pi/theta/2);
    lambda = exp(zz-xx_b);
    df_i = ((lambda+1)./(lambda+c.^2)).^(theta/pi); % df ~ 1/tau
    
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
