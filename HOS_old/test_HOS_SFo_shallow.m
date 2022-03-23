clear
addpath ./HOS_SFo
DO_EXPORT = true;

%% input

ka = 0.1;%.25;
M = 5;
nwt.sim.Nx = 2^9;

lambda = 10;
L = 5*lambda;
h = .15*lambda; 
phaseRad = 0;

k = 2*pi/lambda;

dk = 2*pi/L;
nk0 = k./dk;
assert(mod(nk0,1)==0);


% T = 2*pi/( (1+.5*ka^2)*sqrt( 9.81*k*tanh(k*h) ) );
T = 2*pi/( sqrt( 9.81*k*tanh(k*h) ) );
fprintf('tanh(k*h) = %g.\n',tanh(k*h))

nwt.sim.dt = T;
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







% %%  Plot results
% figure('color','w','Position',[527  0  1056  1000],'name',sprintf('SFo ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime)); hold on;
% nPannel = length(res.t);
% for iP = 1:nPannel
%     subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'k') ;
%     ylabel(sprintf('t = %.2fs',res.t(iP)))
% end

nPannel = size(res.eta,1);
[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814]),[],[0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = zeros(nPannel,1);
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),res.x,res.eta(i,:),'k');
%     ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',res.t(i),nonLinRamp(t_ip(i))))
    ylabel(ha(i),sprintf('t = %.2fs',res.t(i)));
    grid(ha(i),'on');
%     axis(ha(i),'equal')
end
% linkaxes(ha)
xlabel(ha(nPannel),'x [m]','fontsize',11)


if DO_EXPORT
    fileName = sprintf('SFoka%.2gH%.2f',ka,h); fileName(fileName=='.')='p';
    savefig(hf,['./HOS_deepConformalWater/figures/',fileName]);
    export_fig(hf,['./HOS_deepConformalWater/figures/',fileName],'-pdf','-png');
end

%     
% fileName = sprintf('%s%ska%.2g_M%d_H%.2f_Nw%d_dt%.3gT_nx%d',exportPrefix,ka,M,H,NWaves,NT_dt,nx); fileName(fileName=='.')='p';
% if DO_EXPORT
%     copyfile('./proto.m',[exportPath,'/script_',fileName,'.m']) 
%     savefig(hf,[exportPath,'/',fileName]);
%     export_fig(hf,[exportPath,'/',fileName],'-pdf','-png');
% end
% if EXPORT_MAT, save([exportPath,'/',fileName]); end
% 
% 

% %%  Plot results ontop of AHA
% nPannel = length(res.t);
% for iP = 1:nPannel
%     subplot(nPannel,1,nPannel-iP+1), plot(res.x,res.eta(iP,:),'--r') ;
% end

