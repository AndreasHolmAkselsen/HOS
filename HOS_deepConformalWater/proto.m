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
NWaves = 2;
lambda = 10;
ka = .28;
surfaceMethod = 'decayingConformal'; % 'phiS_x', 'phi_x', 'Taylor', 'decayingConformal'
% surfaceMethod = 'Taylor'; % 'phiS_x', 'phi_x', 'Taylor', 'decayingConformal'

L = NWaves*lambda;


% some computations...
g = 9.81;
k0 = 2*pi/lambda;
omega = (1+.5*ka^2)*sqrt(g*k0);
T = 2*pi/omega;
c_p = 2*pi/T/k0;

% Simulation/plotting time
NT_dt = .3;
dt = NT_dt*T;
t_end = 9*dt;

% stability
Tramp = 0;
% Tramp = 1*T;
nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
k_cut = (M+5)*k0;

initialStepODE = 1e-3*T;


%% Plot background current
dx = L/nx;
x0 = (0:nx-1)'*dx;

%% Simulation
xk0 = k0.*x;
phaseAng = 30*pi/180;
% switch initialCondition
%     case 'Stokes3'
% %         eq. 6 and 7 in HOS-memo
%         phiS = ka.*omega/k0^2*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0)));
%         eta = ka/k0*(cos(xk0)+.5*ka*cos(2*xk0)+3/8*ka^2*(cos(3*xk0)-cos(xk0)));
%         % phi0 = ka.*omega/k0^2.*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0))); %.*exp(k0*z)
%         % [~,psi] = getStreamFunction(dx,z,fft(phi0));
%     case 'linearWave'
%         phiS = ka.*sqrt(g*k0)/k0^2*(sin(xk0-phaseAng));
%         eta = ka/k0*(cos(xk0-phaseAng));
%     case 'wavePacket'
%         packet = exp(-(min(abs(x-x0),L-abs(x-x0))/packetWidth).^2);
%         phiS = ka.*sqrt(g*k0)/k0^2*(sin(xk0)).*packet;
%         eta = ka/k0*(cos(xk0)).*packet; 
%     otherwise
%         
% end

phiS0 = ka.*sqrt(g*k0)/k0^2*(sin(xk0-phaseAng));
eta0 = ka/k0*(cos(xk0-phaseAng));
if strcmp(surfaceMethod,'decayingConformal')
    eta_adj = initializeInitCond(x,eta0,5);
    xi = x; kx = getKx(xi);
    FFTeta = fft(eta_adj)/nx;
    f =  xi + 2i*fft(conj(fft(eta_adj)/nx).*(abs(kx)<k_cut&kx>0),[],1)+1i*FFTeta(1);
    phiS_adj = ka.*sqrt(g*k0)/k0^2*(sin(real(f)*k0-phaseAng));
    
    % figure('color','w');
    % f0 =  xi + 2i*fft(conj(fft(eta0)/nx).*(abs(kx)<k_cut&kx>0),[],1);
    % subplot(2,1,1); plot(x,eta0,'-',real(f),imag(f),'--',real(f0),imag(f0),':','linewidth',1.5);ylabel('\eta');title('IC verification')
    % subplot(2,1,2); plot(x,phiS0,'-',real(f),phiS_adj,'--','linewidth',1.5);ylabel('\phi^S');
else
    eta_adj=eta0;    phiS_adj=phiS0;
end

% phiS_ip = phiS;
% eta_ip = eta;


ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);
tic
[t,y] = ode45(@HOSODE45 ,[0,t_end],[phiS_adj;eta_adj],ODEoptions);
% [t,y] = ode45(@HOSODE45 ,0:dt:t_end,[phiS;eta],ODEoptions);
CPUTime = toc;
fprintf('CPU time (AHA): %gs\n',CPUTime);
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



% hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('%s Tramp%g ka=%.3g,M=%d,CPU=%.3g',surfaceMethod,Tramp,ka,M,CPUTime));%[-1587 511 560 1000]
% for iP = 1:nPannel
%     ha(iP) = subplot(nPannel,1,nPannel-iP+1); plot(x_ip(:,iP),eta_ip(:,iP),'k');hold on
%     ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))))
%     box off; grid on;
% end
% linkaxes(ha,'xy');
% set(ha(2:end),'XTick',[]);
% xlim(x([1,end]));



[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814]),[],[0,0]);
ha=flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    plot(ha(i),x_ip(:,i),eta_ip(:,i),'k');
    ylabel(ha(i),sprintf('%.1fs',t_ip(i)));
    grid(ha(i),'on');
end
linkaxes(ha)
% ylim(max(res.eta(:))*[-1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)
    xlim(x([1,end]));


fileName = sprintf('%s%ska%.2g_M%d_Nw%d_dt%.3gT_nx%d',exportPrefix,surfaceMethod,ka,M,NWaves,NT_dt,nx); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./proto.m',[exportPath,'/script_',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf','-png');
end
if EXPORT_MAT, save([exportPath,'/',fileName]); end


