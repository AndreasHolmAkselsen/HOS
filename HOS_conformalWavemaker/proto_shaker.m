clear 
global timeReached 
timeReached = 0; 

g = 9.81;
% g = 1.0;

%% input
h = 5; % Depth
wbl = 1; % hinge depth
wbOverWater = .5; % flap extention obove quiescent waterline
Xmax = 1.01;
N_periods_BM_ramp = 2;
waveRampType = 'none'; % {'linear','tanh','none'}
N_wavelenghtsInL = 10;

% interpolation grid for conformal map
nArrX_far = 200;   % # horizontal points (far field)
nArrX_near = 200;  % # horizontal points (near field)
nArrYDown = 200;    % # vertical points
nTheta = 101;       % # flap angles (odd)

% fz = @fz_yySymmetric;
% fz = @fz_volumeConverving;

% for wave breaking example
param.DO_PADDING = 0;
RK4dt = 0;%2.5e-3; % set to zero to use ODE45
%  and  NT_dt =  5 /9/T; lambda = 2*pi; g=1;


relTolODE = 1e-4;% 1e-8;
N_SSGW = 2^12; % number of modes in SSGW solution

% Plot & export options
DO_EXPORT = 0;
EXPORT_MAT = 0;
PLOT_MAP = 0;
exportPrefix = 'BM_';
exportPath = './figures/';
exportFormatsMap = {'-pdf','-png'};
exportFormats = {'-png','-pdf','-m2'};


%% Simulation
nx = 2^10;

% % specifying period
T = 1.0; omega=2*pi/T; 
kTemp = findWaveNumbers(omega,h(1),0,0);
L = 2*pi/kTemp*N_wavelenghtsInL;


flapPlotTimes = [];%linspace(0,T,50); % leave empty skips flap plotting

% Simulation/plotting time
NT_dt = .25;
dt = NT_dt*T;
param.t_end = 9*dt;
    

% TRamp = 0*T;
% param.nonLinRamp = @(t) max(0,1-exp(-(t/TRamp)^2));
% param.iModeCut = 2*(NWaves + (param.M+5));


param.M = 5; 

% param.iModeCut = inf;
% param.kd__kmax = .5;
% rDamping = .25;

param.iModeCut = inf* (param.M+5)*N_wavelenghtsInL; %nx__wave*NWaves/4;
param.kd__kmax = 0;
rDamping = 0;



dx = L/nx;
% x = (0:nx-1)'*dx;
x = linspace(0,L,nx+1)';x(end)=[];
fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-param.iModeCut/ (nx/2),0) )

t0 = 0;
initialStepODE = 1e-3*T;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);%,'MaxStep',1e-4);

fileName = sprintf('%sT%.2f_M%d_h%.2f_wbl%.2f_Xmax%.1f_L%.3g_dt%.3gT_nx%d_pad%d_ikCut%.4g_Md%.2g_r%.2g',exportPrefix,T,param.M,h,wbl,Xmax,L,NT_dt,nx,param.DO_PADDING,param.iModeCut,param.kd__kmax,rDamping); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./proto_beachRamp.m',[exportPath,'/m/',fileName,'.m']) 
end




%%  map preparation

% create a ramp for the wavemaker motion
Tramp = N_wavelenghtsInL*T;
switch waveRampType
    case 'tanh'
        tanhCutOff = 0.01;
        shapeFactor = -2*atanh(2*tanhCutOff-1); 
        ramp_t0 = .5*shapeFactor*sech(.5*shapeFactor)^2/Tramp;
        t_add = tanhCutOff/ramp_t0;
        ramp = @(t) (t>t_add).*.5.*(1+tanh( shapeFactor*((t-t_add)/Tramp-.5))) + ramp_t0.*t.*(t<=t_add);
        ramp_t = @(t) (t>t_add).*.5.*shapeFactor.*sech(shapeFactor*((t-t_add)/Tramp-.5)).^2/Tramp + ramp_t0.*(t<=t_add);
    case 'linear'
        ramp = @(t) min( t/Tramp, 1);
        ramp_t = @(t) (t<Tramp)/Tramp;
    case 'none'
        ramp = @(t) 1;
        ramp_t = @(t) 0;
    otherwise
        error('waveRampType ''%s'' not recognized.',waveRampType);
end

% xx = x;
% yy = y;

%     zz = xx + 1i*yy;
%     f_zz = 1+0*zz;
%     f = @(zz,theta) zz + Xmax*sin(omega*t);
%     f_t = Xmax*omega*cos(omega*t);


fzIP = @(zz,t) zz + Xmax*sin(omega*t);
map.fy = @(zz,t) imag(zz);
map.fJInv = @(zz,t) 1+0*zz;
map.ft__fzz = @(zz,t) Xmax*omega*cos(omega*t) + 0*zz;


%% test interpolation with plot
if ~isempty(flapPlotTimes)
    xxEnd = 2*h;
    zz_r = linspace(0,xxEnd,30 )+1i*linspace(-h,wbOverWater,200).';
    zz_b = linspace(0,xxEnd,200)+1i*linspace(-h,wbOverWater,10 ).';

    figure('color','w');hold on; grid on
    plot(0,-wbl,'ok','markersize',8,'linewidth',2);
%     axis([-(wbl+wbOverWater)*sin(thetaMax),xxEnd,-h,wbOverWater])
    for iPlot = 1:length(flapPlotTimes), tPlot=flapPlotTimes(iPlot);
        if iPlot>1, delete(hl); end
        hl = [plot(fzIP(zz_r,tPlot),'r','linewidth',1)  
              plot(fzIP(zz_b,tPlot).','b','linewidth',1)
%               plot(fzIP(linspace(0,xxEnd,30 ),tPlot).','.-k','linewidth',1)
          %    plot([0,0,-(wbl+wbOverWater)*tan(thetaMax*sin(omega*tPlot))],[-h,-wbl,wbOverWater]   ,'-k','linewidth',2)
          ];
        drawnow
    end
end


% NB! may need altering
map.xi = x;
map.zzDepth = h;


[varphiS0,eta0_xiReg] = deal(zeros(size(x)));



param.map = map;
param.g = g;
param.rDampingDim = rDamping*2*pi*sqrt(g/L); % check!

%% Run simulation
tic
if RK4dt~=0
    [t,y] = RK4(@(t,Y) HOS_Taylor_closed(t,Y,param) ,[t0,RK4dt,param.t_end],[varphiS0;eta0_xiReg]);
else
    [t,y] = ode45(@(t,Y) HOS_Taylor_closed(t,Y,param) ,[t0,param.t_end],[varphiS0;eta0_xiReg],ODEoptions);
end
fprintf('CPU time: %gs\n',toc);
varphiS = y(:,1:nx); eta = y(:,nx+1:2*nx);

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
% set(ha,'XLim',[minh,maxh],'YLim',[min(imag(zS_ip(:))),max(imag(zS_ip(:)))])
% set(ha,'DataAspectRatio',[1,1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)

% estimate of reflected wave and add it to export file name
% iRefPan = 9;
% ii = real(zS_ip(:,iRefPan))>-.25*L & real(zS_ip(:,iRefPan))<0;
% aRef = .5*( max(imag(zS_ip(ii,iRefPan)))-min(imag(zS_ip(ii,iRefPan))));
% fileName = sprintf('%s_a1%.4f_aRef%.5f',fileName,ka/k1,aRef); fileName(fileName=='.')='p';

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

