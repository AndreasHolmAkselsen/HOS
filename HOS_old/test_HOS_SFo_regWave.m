clear
addpath ./HOS_SFo

%% input
DO_EXPORT = true;
exportPrefix = 'SFo_';
exportPath = 'C:\gits\wave-current\conformalMapping\conformalHOS_current\figures\';


ka = 0.3;
M = 5;
nwt.sim.Nx = 2^8;
L = 100;
lambda = 10;
h = 100;%2*lambda; % 5
phaseRad = 0;

k = 2*pi/lambda;

dk = 2*pi/L;
nk0 = k./dk;
assert(mod(nk0,1)==0);

% T = 2*pi/( (1+.5*ka^2)*sqrt( 9.81*k*tanh(k*h) ) );
T = 2*pi/sqrt(9.81*k*tanh(k*h));
fprintf('tanh(k*h) = %g.\n',tanh(k*h))

dt = 5*T;
KEEP_ALL_TIMES = true;
nwt.sim.tMax = 9*dt;

nwt.solver.LPfilter.type = 'cut'; %'power'
% nwt.solver.LPfilter.kCutMode = nk0*(5+M);
% nwt.solver.LPfilter.kCutMode = 50; % from figure caption
nwt.solver.LPfilter.kCutMode = 1e12; %i.e., inf
nwt.solver.rTol = 1e-4;%1e-8;
nwt.solver.ramp.type = 'exp';
% ramp: F=1-exp(-(t/Tramp)^nRamp);
nwt.solver.ramp.Ta = 1*T; % TRamp
nwt.solver.ramp.n = 2; % nRamp

%% Run
nwt.solver.M = M;
nwt.sim.Lx = L;
nwt.sim.depth = h;
nwt.solver.type = 'hos'; 
nwt.sim.type = 'periodicDomain2D';
wt.init.waveComp{1}.type = 'Regular_ModeNumber';
nwt.init.waveComp{1}.type = 'Regular_ModeNumber';
nwt.init.waveComp{1}.nk0 = nk0;
nwt.init.waveComp{1}.a = ka/k;
nwt.init.waveComp{1}.phaseRad = phaseRad;

if KEEP_ALL_TIMES
    nwt.sim.dt = 0;
else
    nwt.sim.dt = dt;
end

tic
res=runNWT(nwt);
CPUTime = toc;
fprintf('CPU time (SFo): %gs\n',CPUTime);


%%  Plot results
np = floor(nwt.sim.tMax/dt);
if KEEP_ALL_TIMES
    t_ip = linspace(0,res.t(end),np)';
    eta  = interp1(res.t,res.eta ,t_ip);
    t = t_ip; 
    dt = t(end)/(np-1);
else
    eta = res.eta; t = res.t; 
end


[hf, ha] = multi_axes(np,1,figure('color','w','position',[527  0  1056  1000],'name',sprintf('SFo ka=%.3g,M=%d',ka,M)),[],[0,0]);
ha = flipud(ha);
set([ha(2:np).XAxis],'Visible','off');
for iP=1:np
    plot(ha(iP),res.x,eta(iP,:),'k');
    ylabel(ha(iP),sprintf('%.1fs',t(iP)));
    grid(ha(iP),'on');
end
linkaxes(ha)
ylim(max(eta(:))*[-1,1])
xlabel(ha(np),'x [m]','fontsize',11)

if DO_EXPORT
    fileName = sprintf('%ska%.2g_M%d_h%.2f_nx%d_ikCut%d',exportPrefix,ka,M,h,nwt.sim.Nx,nwt.solver.LPfilter.kCutMode); fileName(fileName=='.')='p';
    
    copyfile('./test_HOS_SFo_regWave.m',[exportPath,'/script_',fileName,'.m'])
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf','-png');
end