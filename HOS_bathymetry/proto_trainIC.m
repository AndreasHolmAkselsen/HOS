clear 
global timeReached 
timeReached = 0; 

g = 9.81;
% g = 1.0;

%% input


% interpolation grid for conformal map
map.nArrX_far = 1500;   % # horizontal points (far field)
map.nArrX_near = 500;  % # horizontal points (near field)
map.nArrYDown = 100;    % # vertical points


% Plot & export options
DO_EXPORT = 1;
EXPORT_MAT = 2;
PLOT_MAP = 0;
PLOT_INTERPOLATION_MAP = 0;
exportPrefix =  'train_';
exportPath = './figures/';
exportFormatsMap = {'-pdf','-png'};
exportFormats = {'-png','-pdf','-m2'};



% Domain parameters
boundaryType = 'open'; % 'open' or 'closed'
% NWaves = 5;
% Lx = (2*pi/IC.k)*NWaves;
% nx__wave = 2^8;
% Lx = 75;
% NWaves = 1;
% nx__wave = 1024*1.5;
% xLR = [0,Lx];%[-Lx/2,Lx/2];

% reduced domain
Lx = 75;
N = 1*1024; % N is the number spacings or the number of poits excluding x=+L, i.e., x(end)+dx=L
xLR = [0,Lx];



waveMaker =  []; % no wavemaker

% irregular, slope bathymetry
addpath c:/gits/timsas2/matlabLibs/
map.H = [1,0.5]; % plateau levels
map.xx_b = Lx/2*pi ;  % xi-coordinate of "begining of" edge
map.theta = [90]*pi/180; % slope angles (positive values)
param.nonLinRamp = @(t) 1;
param.iModeCut = N/3;
param.kd__kmax = .0;
rDamping = .0;



% beach
beach.length = 0*Lx/3;
beach.absorption = 1.0;
beach.uniformFraction = 2/3; 

% Initial conditions
IC.depth = map.H(1);
IC.ka = .2; % linear wave steepness
IC.k = pi;%findWaveNumbers(2*pi/IC.T,map.IC.depth,0,0);
IC.T = 2*pi/sqrt(g*IC.k*tanh(IC.k*IC.depth));
INIT_WAVE_TYPE = 'SSGW';  % 'SSGW', 'linear' or a file name
IC.N_SSGW = 2^12; % number of modes in SSGW solution
packetWidth__L = .5;  % set to inf if not simulating wave packets
packetCentre__L = .25;
packetType = 'taperedTrain'; %'taperedTrain', 'gaussian'



% % Simulation/plotting time
% NT_dt = 1;
% dt = NT_dt*IC.T;
% param.t_end = 9*dt;

% Simulation/plotting time
% param.t_end = 'fromBM'; tFinalStill = 20;
NT_dt = 2.5;
param.t_end = 9*NT_dt*IC.T;

% numerical
param.M = 5; 
% initialStepODE = 1e-3*IC.T;
initialStepODE = 1e-3;
relTolODE = 1e-4;% 1e-8;


param.DO_PADDING = 0;
RK4dt = 0;%2.5e-3; % set to zero to use ODE45
%  and  NT_dt =  5 /9/T; lambda = 2*pi; g=1;




%% Simulation
dx = Lx/N;
N_end = N-strcmp(boundaryType,'open');
x = (0:N_end)'*dx+xLR(1);
fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-param.iModeCut/(N/(1+strcmp(boundaryType,'open'))),0) )

switch packetType
    case 'gaussian'
        packet = exp(-((x/Lx-packetCentre__L)/packetWidth__L).^2);
    case 'taperedTrain'
        taperWidth = 2*2*pi/IC.k; % two wavelengths
        packet=zeros(size(x));
        x1 = (packetCentre__L - .5*packetWidth__L)*Lx;
        x2 = (packetCentre__L + .5*packetWidth__L)*Lx;
        ind1=find(x <= x1,1,'last');
        ind2=find(x >= x2,1,'first');
        nRamp = round(taperWidth/dx);
        packet(ind1+(0:nRamp-1)) = .5*(1-cos(linspace(0,pi,nRamp)));
        packet(ind1+nRamp:ind2-nRamp) = 1;
        packet(ind2-(0:nRamp-1))=0.5*(1+cos(linspace(pi,0,nRamp)));
end



Hstr = sprintf('%.2f_',map.H); thetaStr = sprintf('%.0f_',map.theta*180/pi);
% fileName = sprintf('%s%s_%s_T%.2f_ka%.2g_M%d_H%stheta%sNw%d_dt%.3gT_nx%d_pad%d_ikCut%.4g_Md%.2g_r%.2g',exportPrefix,boundaryType,INIT_WAVE_TYPE,IC.T,IC.ka,param.M,Hstr,thetaStr,NWaves,NT_dt,N,param.DO_PADDING,param.iModeCut,param.kd__kmax,rDamping); fileName(fileName=='.')='p';
fileName = sprintf('%s%s_M%d_H%stheta%sNw%.3g_nx%d_L%.3g_Lb%.3g_pad%d_ikCut%.4g_Md%.2g_r%.2g',exportPrefix,boundaryType,param.M,Hstr,thetaStr,NT_dt,N,Lx,beach.length,param.DO_PADDING,param.iModeCut,param.kd__kmax,rDamping); 
if IC.ka>0 
    fileName = sprintf('%s_%s_T%.2f_ka%.2g',fileName,INIT_WAVE_TYPE,IC.T,IC.ka);
end
fileName(fileName=='.')='p';

if DO_EXPORT
    copyfile('./proto_trainIC.m',[exportPath,'/m/',fileName,'.m']) 
end


%% Initial conditions:
switch INIT_WAVE_TYPE
    case 'linear'
        phaseAngs = IC.k.*x-0*pi/180;
        omega = 2*pi/IC.T;
        h0 = IC.ka/IC.k*(cos(phaseAngs));
        phiS0 = IC.ka/IC.k.*g/omega*sin(phaseAngs);
    case 'SSGW'
        [h0,phiS0] = initSSGW(IC.k,IC.depth,IC.ka,IC.N_SSGW,x,g);
    otherwise % input data file assumed
        filePath = ['./IC/',INIT_WAVE_TYPE,'.mat'];
        assert(isfile(filePath));
        load(filePath)
        assert(x0(end)<=x(end))
        h0 = [h0;zeros(N-length(x0),1)];
        varphiS0 = [varphiS0;zeros(N-length(x0),1)];
%         figure, plot(x,eta,x,phiS,'--')
end
% Adjust in case of wave packets.
h0 = h0.*packet;
phiS0 = phiS0.*packet;

% dim.L  = Lx/(2*pi);
% dim.t = sqrt(dim.L/g);
% dim.phi = sqrt(dim.L^3*g);
% dim.U = dim.L/dim.t;


%%  map preparation
map.xLR = xLR;
if ~isempty(map.xx_b)
    [map,fIP,varphiS0,eta0_xiReg,hf_map] = mapDomainSC(map,x,h0,phiS0,PLOT_MAP,PLOT_INTERPOLATION_MAP);
    
    if PLOT_MAP && DO_EXPORT
        Hstr = sprintf('%.2f_',map.H); thetaStr = sprintf('%.0f_',map.theta*180/pi);
        fileNameMap = sprintf('%s%s_ka%.2g_H%stheta%sNw%d',exportPrefix,INIT_WAVE_TYPE,IC.ka,Hstr,thetaStr,NWaves); fileNameMap(fileNameMap=='.')='p';
        export_fig(hf_map(1),['./figures/map/map_',fileNameMap],exportFormatsMap{:})
        savefig(hf_map(1),['./figures/fig/map_',fileNameMap])
        export_fig(hf_map(2),['./figures/map/mapZoom_',fileNameMap],exportFormatsMap{:})
        savefig(hf_map(2),['./figures/fig/mapZoom_',fileNameMap])
    end
else
   assert(isempty(map.theta) && isscalar(map.H));
   map.xi = x;
   map.zzDepth = map.H;
   varphiS0 = phiS0;
   eta0_xiReg = h0;
   map.fy = @(zz) imag(zz);
%    map.fJInv = @(zz) 1;
   map.f_zz = @(zz) ones(size(zz));
   map.xxLR = xLR;
   fIP = @(zz) zz;
   map.zzRoots = nan;
end
x_xiReg = real(fIP(map.xi.')).';

% interpolate as needed for anti-aliasing (padding)
p = 1+3*param.DO_PADDING; Nd=N_end*(p+1)/2;
x_AA = (0:Nd)'*N_end/Nd*dx + x(1); % =linspace(xLR(1),xLR(2),Nd+1)';
x_xiReg_AA = interp1((0:N_end)'/N_end,x_xiReg,(0:Nd)'/Nd);
assert(abs(x_xiReg_AA(1)-x(1))<1e-9 && abs(x_xiReg_AA(end)-x(end))<1e-9)
x_xiReg_AA([1,end]) = x([1,end]);




%% Beach init.
if exist('beach','var') && beach.length>0
    xB = xLR(2)-beach.length;    
    xU = xB + beach.length*(1-beach.uniformFraction);
%     u = (x_xiReg_AA-xB)/(xLR(2)-xB).*(x_xiReg_AA>xB);
    u = (x_xiReg_AA-xB)/(xU-xB).*(x_xiReg_AA>xB);
%     param.beach.nu = beach.absorption*u.^2; % Bonnefoy (2005)/SFo
    param.beach.nu = beach.absorption*u.^2.*(3-2*u); %Bonnefoy (2006)   
    param.beach.nu(x_xiReg_AA>=xU) =  beach.absorption;
%     figure, plot(x_xiReg_AA,param.beach.nu,'.-')
end


%% Linear wavemaker init, copied from SFo git hosm-nwt2d
if ~isempty(waveMaker)
    assert(strcmp(boundaryType,'closed'),'boundaryType must be ''closed'' when specifying a wavemaker load.')
    [phiAdd_of_x,waveMakerTime,X,z] = BM.callSFoWavemakerFunctions(waveMaker,N,Lx,param.DO_PADDING);
    param.waveMaker.time = waveMakerTime;
    param.waveMaker.phiAdd = interp1(x_AA,phiAdd_of_x,x_xiReg_AA);
    
    if strcmp(param.t_end,'fromBM')
        dt_ODE = param.waveMaker.time(2);
        param.t_end = param.waveMaker.time(end)+tFinalStill;
        t_ODE = 0:dt_ODE:param.t_end;
    else
        t_ODE = [0,param.t_end];
    end
else
    t_ODE = [0,param.t_end];
end
% figure, plot(x,squeeze(param.waveMaker.phiAdd(:,1,:)),'k','linewidth',1);hold on;grid on


param.map = map;
param.g = g;
% param.rDampingDim = rDamping*2*pi*sqrt(g/diff(map.xxLR)); % check!
param.rDampingDim = rDamping*sqrt(2*pi*g/diff(map.xxLR)); 
%% Run simulation
switch boundaryType
    case 'open'
        HOS_function = @HOS_Taylor_open;
    case 'closed'     
        HOS_function = @HOS_Taylor_closed;
    otherwise
        error('Domain type ''%s'' not recognized.',boundaryType);
end
tic
if RK4dt~=0
    [t,y] = RK4(@(t,Y) HOS_function(t,Y,param) ,[0,RK4dt,t_ODE(end)],[varphiS0;eta0_xiReg]);
else
    ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);
    [t,y] = ode45(@(t,Y) HOS_function(t,Y,param),t_ODE,[varphiS0;eta0_xiReg],ODEoptions);
end
fprintf('CPU time: %gs\n',toc);
varphiS = y(:,1:end/2); eta = y(:,end/2+1:end);

iNaN = find(isnan(varphiS(:,1)),1,'first');
if ~isempty(iNaN), t(iNaN:end)=[]; varphiS(iNaN:end,:)=[]; eta(iNaN:end,:)=[]; end
clear y



% t_ip = (0:dt:param.t_end)';
t_ip = linspace(0,t(end),10).';
% t_ip = linspace(0,.9*t(end),10).';

nPannel = length(t_ip);
varphiS_ip = interp1(t,varphiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

zS_ip = fIP(map.xi+1i*eta_ip);
maxh = max(real(zS_ip(:)));minh = min(real(zS_ip(:)));

[hf, ha] = multi_axes(nPannel,1,figure('color','w','name',sprintf('Conformal; ka=%.3g',IC.ka)),[.075,.04,.05,.05],[.0,0]);
% [hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('Conformal; ka=%.3g',IC.ka)),[.075,.04,.05,.05],[.0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
zSingular = fIP(map.zzRoots);
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),zS_ip(:,i),'k');
%     ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),param.nonLinRamp(t_ip(i))))
    ylabel(ha(i),sprintf('t = %.2fs',t_ip(i)));
    grid(ha(i),'on');
    plot(ha(i),[1;1].*real(zSingular(:)).',[minh;maxh],'--k'); 
end
% axis(ha,'equal','tight')
set(ha,'XLim',[minh,maxh],'YLim',[min(imag(zS_ip(:))),max(imag(zS_ip(:)))])
% set(ha,'DataAspectRatio',[1,1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)
% set(ha,'XLim',[minh,maxh],'YLim',[-.02,.02])


% % estimate of reflected wave and add it to filename
% iRefPan = 9;
% ii = real(zS_ip(:,iRefPan))>-.25*Lx & real(zS_ip(:,iRefPan))<0;
% aRef = .5*( max(imag(zS_ip(ii,iRefPan)))-min(imag(zS_ip(ii,iRefPan))));
% fileName = sprintf('%s_a1%.4f_aRef%.5f',fileName,IC.ka/IC.k,aRef); fileName(fileName=='.')='p';

if DO_EXPORT
    savefig(hf,[exportPath,'/fig/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],exportFormats{:});
end

wh = whos;
switch EXPORT_MAT 
    case  1
        vars = setdiff({wh.name},{'t','y','varphiS','eta','phiAdd_of_x','param'});
        save([exportPath,'/mat/',fileName],vars{:});
    case 2
        vars = setdiff({wh.name},{'y','varphiS','phiAdd_of_x','param'}); % keep eta
        save([exportPath,'/mat/',fileName],vars{:});
    case 3
        save([exportPath,'/mat/',fileName]); 
end



