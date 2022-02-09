clear
addpath ./HOS_SFo_curr

%% input
% Resolution
nx = 2^10;
M = 5; % solution order
relTolODE = 1e-8;

% % Plot & export options
PLOT_CURRENT = true; 
% PLOT_FINAL_VALOCITY_FIELD = false;
PLOT_WITH_SUBPLOT = false;
DO_EXPORT = true;
EXPORT_MATFILE = true;
exportPrefix = 'SFo_basin';
exportPath = './HOS_SFo_curr/figures/';

KEEP_ALL_TIMES = true;

ka0 = .2; % steepness if not limited by flap angle
TLin = 1;


% Wave specification
h = 100; % 5
g = 9.81;

% Basin
L = 130;
theta_max = 10*pi/180; % Flap angle limit
wbl = 3; % hinge depth


% Current
U_curr = .17;
zeta_j = L/2-3.5i ;% object centre
U_j    =  -(.07+U_curr)*1i ;% object strength-- +1:source, -1:sink, 1i:counter-clockwise vortex, -1i ...
nMirror = 3; % number of times the domain is repeated in x.
% zeta_j = [7-3.5i,  50-100i ];% object centre
% U_j    = [  -.25i,   .07i   ];% object strength-- +1:source, -1:sink, 1i:counter-clockwise vortex, -1i ...

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

phaseRad = 0;




% Simulation/plotting time
NT_dt = 15;%5 % Number of periods between plot panel.
nPannel = 10; % Number of time panel to display.
NT_ramp = 1; % Number of periods in nonlinearity ramp.
% Simulation/plotting time
dt = NT_dt*T;
t_end = (nPannel-1)*dt;


% stability
TRamp = 1*T;
rampExponent = 2;
dk = 2*pi/L;
nk0 = k0./dk;
kCutMode = nk0*(5+M);
initialStepODE = 1e-3*T;
nwt.solver.deAliasingProductOrder = 4; % 4

fprintf('tanh(k*h) = %g.\n',tanh(k0*h))


%% Run



exportName = sprintf('%sT%.1f_ka%.2g_M%d_Nw%d_dt%.3gT_L%.0f',exportPrefix,T,ka,M,NWaves,NT_dt,L); exportName(exportName=='.')='p';
if PLOT_CURRENT && ~isempty(U_j)
    dx = L/nx;
    x = (0:nx-1)'*dx;
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
    
    clear x
    if DO_EXPORT
        savefig(hf_c,[exportPath,'/curr_',exportName]);
        export_fig(hf_c,[exportPath,'/curr_',exportName],'-pdf');
        savefig(hf_Usurf,[exportPath,'/Usruf_',exportName]);
        export_fig(hf_Usurf,[exportPath,'/Usruf_',exportName],'-pdf');
    end
    
end


if KEEP_ALL_TIMES
    nwt.sim.dt = 0;
else
    nwt.sim.dt = dt;
end
nwt.sim.tMax = t_end;
nwt.sim.Nx = nx; 
nwt.solver.M = M;
nwt.sim.Lx = L;
nwt.sim.depth = h;
nwt.solver.type = 'hos'; 
nwt.sim.type = 'periodicDomain2D';
wt.init.waveComp{1}.type = 'Regular_ModeNumber';
nwt.init.waveComp{1}.type = 'Regular_ModeNumber';
nwt.init.waveComp{1}.nk0 = nk0;
nwt.init.waveComp{1}.a = ka/k0;
nwt.init.waveComp{1}.phaseRad = phaseRad;


nwt.solver.ramp.type = 'exp';
% ramp: F=1-exp(-(t/Tramp)^nRamp);
nwt.solver.ramp.Ta = TRamp; % TRamp
nwt.solver.ramp.n = rampExponent; % nRamp

nwt.solver.LPfilter.type = 'cut'; %'power'
nwt.solver.LPfilter.kCutMode = kCutMode;
nwt.solver.rTol = relTolODE;
nwt.solver.initialStepODE = initialStepODE;

nwt.curr.df = df;

tic
res=runNWT(nwt);
CPUTime = toc;
fprintf('CPU time (SFo): %gs\n',CPUTime);


% KEEP_ALL_TIMES = false;
% t_ip = linspace(0,res.t(end),nPannel)';
% res.eta  = interp1(res.t,res.eta ,t_ip);
% res.t = t_ip;

%%  Plot results
if KEEP_ALL_TIMES
    % single axes plot
    hf = figure('color','w','Position',[527  0  2056  1500],'name',sprintf('SFo ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
    hold on;box off;
    set(gca,'XAxisLocation','origin');
    a__height = 300;
    dz__a = 1;
    nz = a__height/dz__a;
    zz = dz__a*(0:nz-1)'*max(abs(res.eta(1,:)));
    t_ip = linspace(0,res.t(end),nz)';
    eta_ip  = interp1(res.t,res.eta ,t_ip);
    plot(res.x,eta_ip+zz,'k')
    xlabel('x [m]'); ylabel('\eta [m]');
elseif PLOT_WITH_SUBPLOT
    % subplot
    hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('SFo ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
    nPannel = length(res.t);
    for iP = 1:nPannel
        subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'k');hold on
        ylabel(sprintf('t = %.2fs',res.t(iP)))
        set(gca,'XTick',[],'XAxisLocation','origin');
        box off; grid on;
        % plot marker
        xDrift = mod(res.t(iP)*c_p,L);
        yMarker =  interp1(res.x,res.eta(iP,:),xDrift);
        hold on; plot(xDrift,yMarker,'ok','markersize',10)
        text(xDrift,yMarker,num2str(floor(res.t(iP)*c_p/L)),'fontsize',10,'VerticalAlignment','middle','HorizontalAlignment','center')
    end
else
    % single axes plot
    hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('SFo ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
    hold on;box off;
    set(gca,'XAxisLocation','origin');
    dz = 2.5*max(abs(res.eta(1,:)));
    nPannel = length(res.t);
    for iP = 1:nPannel
        z0 = (iP-1)*dz;
        plot(res.x,res.eta(iP,:)+z0,'k','linewidth',1);
        text(L,z0,sprintf('t = %.2fs',res.t(iP)));
        plot([0,L],[z0,z0],'k','linewidth',.5)
    end
    xDrift = mod(res.t*c_p,L);
    iBreak = find( diff(xDrift)<0 )+1; iBreak = [1;iBreak;length(res.t)];
    for iB=1:length(iBreak)-1, ii = iBreak([iB;iB+1])-[0;1];  plot(xDrift(ii),res.t(ii)/dt,'--r');end
    ylim([min(res.eta(1,:)),2*max(res.eta(1,:))+z0])
    xlabel('x [m]'); ylabel('\eta [m]');
end
drawnow


% %%  Plot results ontop of AHA
% nPannel = length(res.t);
% for iP = 1:nPannel
%     subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'--r') ;
% end

fileName = sprintf('%ska%.2g_M%d_Nw%d_dt%.3gT_nx%d',exportPrefix,ka,M,NWaves,NT_dt,nx); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./test_HOS_SFo_basin.m',[exportPath,'/script_',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf');
end
if EXPORT_MATFILE
    clear hf hf_c
    save([exportPath,'/',fileName])
end