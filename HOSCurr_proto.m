clear
global M h x k_cut nonLinRamp df

%% input
% Resolution
nx = 2^9;
M = 5; % solution order
relTolODE = 1e-8;

% Plot & export options
PLOT_CURRENT = true; 
PLOT_FINAL_VALOCITY_FIELD = true;
DO_EXPORT = false;
EXPORT_MATFILE = false;
exportPrefix = 'AHA_';
exportPath = './HOS_SFo_curr/figures/';

% Wave specification
NWaves = 10;
lambda = 10;
ka = .2;

L = NWaves*lambda;
initialCondition = 'linearWave'; % {'linearWave','Stokes3','wavePacket'} 
% if initialCondition='wavePacket':
% packetWidth = .1*L;
% x0 = 2/5*L;

% some computations...
g = 9.81;
k0 = 2*pi/lambda;
omega = (1+.5*ka^2)*sqrt(g*k0);
T = 2*pi/omega;
c_p = 2*pi/T/k0;

% Simulation/plotting time
NT_dt = 10;
dt = NT_dt*T;
t_end = 9*dt;

% stability
Tramp = 1*T;
nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
k_cut = (M+5)*k0;

initialStepODE = 1e-3*T;

%% Specify background current
zeta_j = []; F_j = [];
nMirror = 3; % number of times the domain is repeated in x.

% % similar to basin, source + vortex
% zeta_j = [-.05-.3i,  .1-.1i ]*L;% object centre
% F_j    = [       6,      2i   ];% object strength-- +1:source, -1:sink, 1i:counter-clockwise vortex, -1i ...


% % % similar to basin, two vortices
% zeta_j = [.6-.7i,  .1-.075i ]*L;% object centre
% F_j    = [  .02i,   -.14i   ];% object strength-- +1:source, -1:sink, 1i:counter-clockwise vortex, -1i ...

% single vortex
zeta_j = [.5-.075i  ]*L;% object centre
F_j    = [ -.2i  ];% object strength

% % source + vortex + sink
% zeta_j = [.25-.1i,.5-.1i ,.75-.1i   ]*L;% object centre
% F_j    = [  .2, -.2i,-.2  ];% object strength

% % single doublet
% zeta_j = [.5-.15i  ]*L;% object centre
% F_j    = [ .1i  ];% object strength


F_j = shiftdim(F_j,-1); zeta_j = shiftdim(zeta_j,-1);% ID_j = shiftdim(ID_j,-1);
zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);


% % vortex/source/sink
A_j = .5*F_j.*c_p.*abs(imag(zeta_j));
f  = @(zeta) sum(A_j.*log(zeta-zeta_j) + conj(A_j.*log(conj(zeta)-zeta_j)),3:4);
df = @(zeta) sum(A_j./(zeta-zeta_j) + conj(A_j./(conj(zeta)-zeta_j)),3:4);

% doublet
% A_j = .5*F_j.*c_p.*abs(imag(zeta_j)).^2;
% f = @(zeta) sum(-A_j./(zeta-zeta_j) - conj(A_j)./(zeta-conj(zeta_j)),3:4);
% df = @(zeta) sum(A_j./(zeta-zeta_j).^2 + conj(A_j)./(zeta-conj(zeta_j)).^2,3:4);


%% Plot background current
dx = L/nx;
x = (0:nx-1)'*dx;

if isempty(F_j), df = @(zeta) 0; end

if PLOT_CURRENT && ~isempty(F_j)
%     z = linspace(-.3*L,.3*L,100);
    z = 0:dx:.3*L; z = [-z(end:-1:2),z];
    hf_c = figure('color','w'); ha = gca;
    % plot velocity intensity |U|
    absU = abs(df(x+1i*z));
    ULim = 2.5*max(abs(df(x)));
    hIm = imagesc(x',z',absU',[0,ULim]); 
%     colorbar
    hIm.AlphaData = (absU<ULim)';
    % plot streamlines
    hold on, axis equal xy
    contour(x',z',imag(f(x+1i*z))',20,'k');
    plot(ha.XLim,[0,0],'k')
    ylim([z(1),0])
    
    % add quiver plot
    zeta_ip = linspace(x(1),x(end),10) + 1i*linspace(z(1),z(end),10)';
    df_ip = df(zeta_ip);
    df_ip(abs(df_ip)>ULim) = nan;
    quiver(real(zeta_ip),imag(zeta_ip),real(df_ip),-imag(df_ip),'r');
    
    ha.Visible='off';
    
%     % plot horizontal velocity at z=0
%     figure('color','w');
%     plot(x,real(df(x)))
%     title(sprintf('|\\phi_x^{(1)}| = %.2gm/s, c_p = %.2gm/s',ka*sqrt(g/k0),c_p));
    
    fprintf('c_p = %.3gm/s, max |U(0)| = %.3gm/s, fraction: %.3g\n',c_p,max(abs(df(x))),max(abs(df(x)))/c_p)
    drawnow
    
    if DO_EXPORT
        fileName = sprintf('curr_%ska%.2g_M%d_Nw%d_dt%.3gT',exportPrefix,ka,M,NWaves,NT_dt); fileName(fileName=='.')='p';
        savefig(hf_c,[exportPath,'/',fileName]);
        export_fig(hf_c,[exportPath,'/',fileName],'-pdf');
    end
    
end
% return

%% Simulation
xk0 = k0.*x;
switch initialCondition
    case 'Stokes3'
%         eq. 6 and 7 in HOS-memo
        phiS = ka.*omega/k0^2*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0)));
        eta = ka/k0*(cos(xk0)+.5*ka*cos(2*xk0)+3/8*ka^2*(cos(3*xk0)-cos(xk0)));
        % phi0 = ka.*omega/k0^2.*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0))); %.*exp(k0*z)
        % [~,psi] = getStreamFunction(dx,z,fft(phi0));
    case 'linearWave'
        phiS = ka.*sqrt(g*k0)/k0^2*(sin(xk0));
        eta = ka/k0*(cos(xk0));
    case 'wavePacket'
        packet = exp(-(min(abs(x-x0),L-abs(x-x0))/packetWidth).^2);
        phiS = ka.*sqrt(g*k0)/k0^2*(sin(xk0)).*packet;
        eta = ka/k0*(cos(xk0)).*packet; 
    otherwise
        
end

ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);
tic
[t,y] = ode45(@HOSODE45 ,[0,t_end],[phiS;eta],ODEoptions);
% [t,y] = ode45(@HOSODE45 ,0:dt:t_end,[phiS;eta],ODEoptions);
CPUTime = toc;
fprintf('CPU time (AHA): %gs\n',CPUTime);
phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);
% interpolate to perscribed times

% t_ip = (0:dt:t_end)';
t_ip = linspace(0,t(end),10)';
nPannel = length(t_ip);
phiS_ip = interp1(t,phiS,t_ip);
eta_ip  = interp1(t,eta ,t_ip);

% nPannel = size(eta,1);
% eta_ip = eta; t_ip = t;

hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('AHA ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
for iP = 1:nPannel
    subplot(nPannel,1,nPannel-iP+1), plot(x,eta_ip(iP,:),'k');hold on
    ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))))
    set(gca,'XTick',[]);
    box off; grid on;
end

fileName = sprintf('%ska%.2g_M%d_Nw%d_dt%.3gT_nx%d',exportPrefix,ka,M,NWaves,NT_dt,nx); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./HOSCurr_proto.m',[exportPath,'/script_',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf');
end
if EXPORT_MATFILE
    clear hf hf_c
    save([exportPath,'/',fileName])
end

if PLOT_FINAL_VALOCITY_FIELD
    phiSEnd = phiS(end,:)'; etaEnd = eta(end,:)';
    z = linspace(10*min(etaEnd),max(etaEnd),500);
    [~,~,~,~,hphi,kx] = phiComponentsHOS(phiSEnd,etaEnd);
    [phi,psi] = getStreamFunction(dx,z,hphi);
    
    figure('color','w');
    contour(x',z',(psi+imag(f(x+1i*z)))')
    hold on; 
    plot(x,etaEnd,'r','linewidth',1);
    
end
