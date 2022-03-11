clear
global M x k_cut nonLinRamp  surfaceMethod timeReached t_end
timeReached = 0;

%% input
% Resolution
nx = 2^9;
M = 5; % solution order
relTolODE = 1e-8;

% Plot & export options
DO_EXPORT = false;
EXPORT_MAT = false;
exportPrefix = '';
exportPath = './figures/';

% Wave specification
NWaves = 6;
lambda = 10;
ka = .25;
% surfaceMethod = 'decayingConformal'; % 'phiS_x', 'phi_x', 'Taylor', 'decayingConformal'
% surfaceMethod = 'Taylor'; % 'phiS_x', 'phi_x', 'Taylor', 'decayingConformal'


L = NWaves*lambda;
g = 9.81;
k0 = 2*pi/lambda;



% some computations...
omega = (1+.5*ka^2)*sqrt(g*k0);
T = 2*pi/omega;
c_p = 2*pi/T/k0;


% Simulation/plotting time
NT_dt = 1.6;
dt = NT_dt*T;
t_end = 9*dt;

    
% stability
Tramp = 1*T;
nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
k_cutTaylor = (M+5)*k0;
k_cut_conformal = (nx*pi/L)/4;

T_init = 2*Tramp;
t0 = 0;


initialStepODE = 1e-3*T;


%% Plot background current
dx = L/nx;
x = (0:nx-1)'*dx;

%% Simulation
xk0 = k0.*x;
phaseAng = 30*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);

phiS0 = ka.*sqrt(g*k0)/k0^2*(sin(xk0-phaseAng));
eta0 = ka/k0*(cos(xk0-phaseAng));
k_cut = k_cutTaylor;

if strcmp(surfaceMethod,'decayingConformal')
    surfaceMethod = 'Taylor';
    tic
    [tInit,yInit] = ode45(@HOSODE45 ,[t0,t0+T_init],[phiS0;eta0],ODEoptions);
    fprintf('CPU time init stage: %gs\n',toc);
    surfaceMethod = 'decayingConformal';
    phiS = yInit(end,1:nx).'; eta = yInit(end,nx+1:2*nx).';
    t0 = tInit(end);
    
    eta_adj = initializeInitCond(x,eta,5);
    xi = x; kx = getKx(xi);
    k_cut = k_cut_conformal;
    FFTeta = fft(eta_adj)/nx;
    f =  xi + 2i*fft(conj(fft(eta_adj)/nx).*(abs(kx)<k_cut&kx>0),[],1)+1i*FFTeta(1);
    phiS_adj = interp1([x-L;x;x+L],[phiS;phiS;phiS],real(f));
%     figure, plot(x,phiS,'-',real(f),phiS_adj,'--','linewidth',1.5)
    
    % figure('color','w');
    % f0 =  xi + 2i*fft(conj(fft(eta0)/nx).*(abs(kx)<k_cut&kx>0),[],1);
    % subplot(2,1,1); plot(x,eta0,'-',real(f),imag(f),'--',real(f0),imag(f0),':','linewidth',1.5);ylabel('\eta');title('IC verification')
    % subplot(2,1,2); plot(x,phiS0,'-',real(f),phiS_adj,'--','linewidth',1.5);ylabel('\phi^S');
else
    [phiS_adj,eta_adj] = deal(phiS0,eta0);
    [tInit,yInit] = deal([]);
end

% phiS_ip = phiS;
% eta_ip = eta;


tic
[t2,y2] = ode45(@HOSODE45 ,[t0,t_end],[phiS_adj;eta_adj],ODEoptions);
fprintf('CPU time: %gs\n',toc);

iNaN = find(isnan(y2(:,1)),1,'first');
if ~isempty(iNaN), t2(iNaN:end)=[]; y2(iNaN:end,:)=[]; end
t = [tInit(1:end-1);t2]; y = [yInit(1:end-1,:);y2];
phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);
% interpolate to perscribed times


% t_ip = (0:dt:t_end)';
t_ip = linspace(0,t(end),10).';
nPannel = length(t_ip);
phiS_ip = interp1(t,phiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

if strcmp(surfaceMethod,'decayingConformal')
    xi = x;
    kx = getKx(x);
    f =  xi + 2i*fft(conj(fft(eta_ip,[],1)/nx).*(abs(kx)<k_cut&kx>0),[],1); 
%     max(abs(imag(f)-eta_ip)) % test 
    
%     it = 4;
%     hf = figure('color','w'); hold on
%     [xi2,sig2]= ndgrid(xi,linspace(-(x(end)/NWaves),0,100));
%     f2 = xi2+1i*sig2 + 2i*fft( conj(fft(eta_ip(:,it),[],1)/nx).*(kx>0).*exp(kx.*sig2) ,[],1); 
%     contour(real(f2),imag(f2),xi2,'r'); hold on
%     contour(real(f2),imag(f2),sig2,'b');
%     plot(real(f(:,it)),imag(f(:,it)),'k','linewidth',2)
%     axis equal

    x_ip = real(f);
%     eta_ip = imag(f); %per definition
else
    x_ip = x+0*eta_ip;
end


[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('%s Tramp%g ka=%.3g,M=%d',surfaceMethod,Tramp,ka,M)),[],[0,0]);
ha=flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    plot(ha(i),x_ip(:,i),eta_ip(:,i),'k');
    ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),nonLinRamp(t_ip(i))))
    grid(ha(i),'on');
%     axis(ha(i),'equal')
end
linkaxes(ha)
% ylim(max(res.eta(:))*[-1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)
xlim(x([1,end]));%set(ha,'XLim',x([1,end]));
    
fileName = sprintf('%s%ska%.2g_M%d_Nw%d_dt%.3gT_nx%d',exportPrefix,surfaceMethod,ka,M,NWaves,NT_dt,nx); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./proto.m',[exportPath,'/script_',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf','-png');
end
if EXPORT_MAT, save([exportPath,'/',fileName]); end




% hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('%s Tramp%g ka=%.3g,M=%d,CPU=%.3g',surfaceMethod,Tramp,ka,M,CPUTime));%[-1587 511 560 1000]
% for iP = 1:nPannel
%     ha(iP) = subplot(nPannel,1,nPannel-iP+1); plot(x_ip(:,iP),eta_ip(:,iP),'k');hold on
%     ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))))
%     box off; grid on;
% end
% linkaxes(ha,'xy');
% set(ha(2:end),'XTick',[]);
% xlim(x([1,end]));