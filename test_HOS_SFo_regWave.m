clear
addpath ./HOS_SFo

%% input

ka = 0.25;
M = 5;
nwt.sim.Nx = 2^10;
L = 100;
lambda = 10;
h = 100; % 5
phaseRad = 0;

k = 2*pi/lambda;

dk = 2*pi/L;
nk0 = k./dk;
assert(mod(nk0,1)==0);

T = 2*pi/( (1+.5*ka^2)*sqrt( 9.81*k*tanh(k*h) ) );
fprintf('tanh(k*h) = %g.\n',tanh(k*h))

nwt.sim.dt = 5*T;
nwt.sim.tMax = 9*nwt.sim.dt;

nwt.solver.LPfilter.type = 'cut'; %'power'
nwt.solver.LPfilter.kCutMode = nk0*(5+M);
nwt.solver.rTol = 1e-8;
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

tic
res=runNWT(nwt);
CPUTime = toc;
fprintf('CPU time (SFo): %gs\n',CPUTime);

%%  Plot results
figure('color','w','Position',[527  0  1056  1000],'name',sprintf('SFo ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime)); hold on;
nPannel = length(res.t);
for iP = 1:nPannel
    subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'k') ;
    ylabel(sprintf('t = %.2fs',res.t(iP)))
end


% %%  Plot results ontop of AHA
% nPannel = length(res.t);
% for iP = 1:nPannel
%     subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'--r') ;
% end

