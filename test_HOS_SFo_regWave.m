clear
addpath ./HOS_SFo

%% input

ka = 0.28;
nwt.sim.depth = 5;
nwt.sim.Nx =  2^10;
nwt.sim.Lx = 100;
nwt.solver.M = 5;

lambda = 10;
k = 2*pi/lambda;

phaseRad = 0;

dk = 2*pi/L;
nk0 = k./dk;
assert(mod(nk0,1)==0);

T = 2*pi/( (1+.5*ka^2)*sqrt( 9.81*k*tanh(k*h) ) );

nwt.sim.dt = 2*T / 5;
nwt.sim.tMax = 9*nwt.sim.dt;

nwt.solver.LPfilter.type = 'cut'; %'power'
nwt.solver.LPfilter.kCutMode = nk0*(5+M);
nwt.solver.rTol = 1e-8;
nwt.solver.ramp.type = 'exp';
% ramp: F=1-exp(-(t/Tramp)^nRamp);
nwt.solver.ramp.Ta = 2*T; % TRamp
nwt.solver.ramp.n = 2; % nRamp

%% Run


nwt.solver.type = 'hos'; 
nwt.sim.type = 'periodicDomain2D';
wt.init.waveComp{1}.type = 'Regular_ModeNumber';
nwt.init.waveComp{1}.type = 'Regular_ModeNumber';
nwt.init.waveComp{1}.nk0 = nk0;
nwt.init.waveComp{1}.a = ka/k;
nwt.init.waveComp{1}.phaseRad = phaseRad;


res=runNWT(nwt);


%%  Plot results
figure('color','w','Position',[-1587 511 560 1000]); hold on;
nPannel = length(res.t);
for iP = 1:nPannel
    subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'k') ;
    ylabel(sprintf('t = %.2fs',res.t(iP)))
end
