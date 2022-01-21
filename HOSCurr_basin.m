clear
global M x k_cut nonLinRamp df
g = 9.81;

%% input
% Resolution
nx = 2^11;
M = 5; % solution order
relTolODE = 1e-8;

% Plot & export options
PLOT_CURRENT = true; 
DO_EXPORT = false;
EXPORT_MATFILE = true;
PLOT_FINAL_VALOCITY_FIELD = false;
PLOT_WITH_SUBPLOT = false;

exportPrefix = '';
exportPath = './doc/figures/basin_L130/';

% Wave
ka0 = .2; % steepness if not limited by flap angle
% TLin = 1.5; %  target peroiod
for TLin = 1:.5:5

% Basin
L = 130;
theta_max = 10*pi/180; % Flap angle limit
wbl = 3; % hinge depth

% Time...
NT_dt = 5; % Number of periods between plot panel.
Npannels = 10; % Number of time panel to display.
NT_ramp = 1; % Number of periods in nonlinearity ramp.

% Current
zeta_j = []; U_j = []; U_curr = 0;
% U_curr = .17;
zeta_j = L/2-3.5i ;% object centre
U_j    =  -(.07+U_curr)*1i ;% object strength-- +1:source, -1:sink, 1i:counter-clockwise vortex, -1i ...
nMirror = 3; % number of times the domain is repeated in x.
% zeta_j = [7-3.5i,  50-100i ];% object centre
% U_j    = [  -.25i,   .07i   ];% object strength-- +1:source, -1:sink, 1i:counter-clockwise vortex, -1i ...





%% Code..
U_j = shiftdim(U_j,-1); zeta_j = shiftdim(zeta_j,-1);% ID_j = shiftdim(ID_j,-1);
zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);


% % vortex/source/sink
A_j = .5*U_j.*abs(imag(zeta_j));
f  = @(zeta) sum(A_j.*log(zeta-zeta_j) + conj(A_j.*log(conj(zeta)-zeta_j)),3:4) + U_curr.*zeta;
df = @(zeta) sum(A_j./(zeta-zeta_j) + conj(A_j./(conj(zeta)-zeta_j)),3:4) + U_curr ;


% some computations...
wLin = 2*pi/TLin;
if U_curr == 0
    kLin = wLin^2/g;
else
    kLin = ( g+2*U_curr*wLin - sqrt(g^2+4*g*U_curr*wLin) )/(2*U_curr^2);
end
NWaves = round(L./(2*pi/kLin));
lambda = L/NWaves;
k0 = 2*pi/lambda;


% limit ka according to a flap angle limit theta_max
c1 = 2-4*sinh(.5*k0*wbl).*exp(-.5*k0*wbl)./(k0*wbl); % deep water limit
ka = min(theta_max*wbl*k0.*c1,ka0);
fprintf('ka used: %g.\n',ka);

T = 2*pi/((1+.5*ka^2)*sqrt(g*k0) + k0*U_curr);
c_p = 2*pi/T/k0;


% Simulation/plotting time
dt = NT_dt*T;
t_end = (Npannels-1)*dt;

% stability
Tramp = NT_ramp*T;
nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
k_cut = (M+5)*k0;
initialStepODE = 1e-3*T;

%% Plot background current
dx = L/nx;
x = (0:nx-1)'*dx;

if isempty(U_j), df = @(zeta) 0; end

exportName = sprintf('%sT%.1f_ka%.2g_M%d_Nw%d_dt%.3gT_L%.0f',exportPrefix,T,ka,M,NWaves,NT_dt,L); exportName(exportName=='.')='p';
if PLOT_CURRENT && ~isempty(U_j)
    z = 0:dx:15; z = [-z(end:-1:2),z];
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
%     xlim(xLim)
    
%     % add quiver plot
%     zeta_ip = linspace(x(1),x(end),10) + 1i*linspace(z(1),z(end),10)';
%     df_ip = df(zeta_ip);
%     df_ip(abs(df_ip)>ULim) = nan;
%     quiver(real(zeta_ip),imag(zeta_ip),real(df_ip),-imag(df_ip),'r');
    
    ha.Visible='off';
    
    % plot horizontal velocity at z=0
    hf_Usurf = figure('color','w','position',[2369 747 560 178]);
    plot(x,real(df(x)))
    title(sprintf('|\\phi_x^{(1)}| = %.2gm/s, c_p = %.2gm/s',ka*sqrt(g/k0),c_p));
    ylabel('\Phi_x(x,0)');xlabel('x'); grid on;
    xlim(x([1,end]))
    
    fprintf('c_p = %.3gm/s, max |U(0)| = %.3gm/s, fraction: %.3g\n',c_p,max(abs(df(x))),max(abs(df(x)))/c_p)
    drawnow
    
    if DO_EXPORT
        savefig(hf_c,[exportPath,'/curr_',exportName]);
        export_fig(hf_c,[exportPath,'/curr_',exportName],'-pdf');
        savefig(hf_Usurf,[exportPath,'/Usruf_',exportName]);
        export_fig(hf_Usurf,[exportPath,'/Usruf_',exportName],'-pdf');
    end
    
end

%% Simulation
xk0 = k0.*x;

initialCondition = 'linearWave'; % {'linearWave','Stokes3','wavePacket'} 

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
CPUTime = toc;
fprintf('CPU time (AHA): %gs\n',CPUTime);
phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);

% interpolate to perscribed times
% t_ip = linspace(0,t(end),Npannels)';
dt = t(end)/(Npannels-1);
t_ip = (0:dt:t(end))';
nPannel = length(t_ip);
% phiS_ip = interp1(t,phiS,t_ip);
eta_ip  = interp1(t,eta ,t_ip);
hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('AHA ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
if PLOT_WITH_SUBPLOT
    % subplot
    for iP = 1:nPannel
        subplot(nPannel,1,nPannel-iP+1), plot(x,eta_ip(iP,:),'k');hold on
        ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))))
        set(gca,'XTick',[],'XAxisLocation','origin');
        box off; grid on;
        % plot marker
        xDrift = mod(t_ip(iP)*c_p,L);
        yMarker =  interp1(x,eta_ip(iP,:),xDrift);
        hold on; plot(xDrift,yMarker,'ok','markersize',10)
        text(xDrift,yMarker,num2str(floor(t_ip(iP)*c_p/L)),'fontsize',10,'VerticalAlignment','middle','HorizontalAlignment','center')
    end
else
    % single axes plot
    hold on;box off;
    set(gca,'XAxisLocation','origin');
    dz = 2.5*max(abs(eta_ip(1,:)));
    for iP = 1:nPannel
        z0 = (iP-1)*dz;
        plot(x,eta_ip(iP,:)+z0,'k','linewidth',1);
        text(L,z0,sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))));
        plot([0,L],[z0,z0],'k','linewidth',.5)
    end
    xDrift = mod(t*c_p,L);
    iBreak = find( diff(xDrift)<0 )+1; iBreak = [1;iBreak;length(t)];
    for iB=1:length(iBreak)-1, ii = iBreak([iB;iB+1])-[0;1];  plot(xDrift(ii),t(ii)/dt,'--r');end
    ylim([min(eta_ip(1,:)),2*max(eta_ip(1,:))+z0])
    xlabel('x [m]'); ylabel('\eta [m]');
end

if DO_EXPORT
    copyfile('./HOSCurr_proto.m',[exportPath,'/script_',exportName,'.m']) 
    savefig(hf,[exportPath,'/',exportName]);
    export_fig(hf,[exportPath,'/',exportName],'-pdf');
end
if EXPORT_MATFILE
    clear hf hf_c
    save([exportPath,'/',exportName])
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

end