clear
addpath ./HOS_SFo_curr

%% input
% Resolution
nx = 2^8;
M = 5; % solution order
relTolODE = 1e-8;

% % Plot & export options
% PLOT_CURRENT = false; 
% PLOT_FINAL_VALOCITY_FIELD = false;
DO_EXPORT = false;
EXPORT_MATFILE = false;
exportPrefix = 'SFo_';
exportPath = './HOS_SFo_curr/figures/';

% Wave specification
NWaves = 10;
lambda = 10;
ka = .28;
h = 100; % 5

L = NWaves*lambda;


% some computations...
g = 9.81;
k0 = 2*pi/lambda;
omega = (1+.5*ka^2)*sqrt(g*k0);
T = 2*pi/omega;
c_p = 2*pi/T/k0;
dk = 2*pi/L;
nk0 = k0./dk;

phaseRad = 0;


% Simulation/plotting time
NT_dt = 10;
dt = NT_dt*T;
t_end = 9*dt;

% NT_dt=0;dt = 0;


% curr
zeta_j = [.5-.075i  ]*L;% object centre
F_j    = 0*[ -.2i  ];% object strength
nMirror = 3; % number of times the domain is repeated in x.


% stability
TRamp = 1*T;
rampExponent = 2;
kCutMode = nk0*(5+M);
initialStepODE = 1e-3*T;
nwt.solver.deAliasingProductOrder = 4; % 4

fprintf('tanh(k*h) = %g.\n',tanh(k0*h))




%% Run
nwt.sim.dt = dt;
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

% single vortex
F_j = shiftdim(F_j,-1); zeta_j = shiftdim(zeta_j,-1);% ID_j = shiftdim(ID_j,-1);
zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);

% % vortex/source/sink
A_j = .5*F_j.*c_p.*abs(imag(zeta_j));
f  = @(zeta) sum(A_j.*log(zeta-zeta_j) + conj(A_j.*log(conj(zeta)-zeta_j)),3:4);
df = @(zeta) sum(A_j./(zeta-zeta_j) + conj(A_j./(conj(zeta)-zeta_j)),3:4);

nwt.curr.df = df;


tic
res=runNWT(nwt);
CPUTime = toc;
fprintf('CPU time (SFo): %gs\n',CPUTime);

%%  Plot results
hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('SFo ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime)); hold on;
if dt == 0
    subplot(3,1,1), plot(res.x,res.eta(round(end/4),:),'k') ;
    ylabel(sprintf('t = %.2fs',res.t(round(end/4))))
    subplot(3,1,2), plot(res.x,res.eta(round(end/2),:),'k') ;
    ylabel(sprintf('t = %.2fs',res.t(round(end/2))))
    subplot(3,1,3), plot(res.x,res.eta(end,:),'k') ;
    ylabel(sprintf('t = %.2fs',res.t(end)))
else
    nPannel = length(res.t);
    for iP = 1:nPannel
        subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'k') ;
        ylabel(sprintf('t = %.2fs',res.t(iP)))
    end
end


% %%  Plot results ontop of AHA
% nPannel = length(res.t);
% for iP = 1:nPannel
%     subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'--r') ;
% end

fileName = sprintf('%ska%.2g_M%d_Nw%d_dt%.3gT_nx%d',exportPrefix,ka,M,NWaves,NT_dt,nx); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./test_HOS_SFo_Curr_regWave.m',[exportPath,'/script_',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf');
end
if EXPORT_MATFILE
    clear hf hf_c
    save([exportPath,'/',fileName])
end