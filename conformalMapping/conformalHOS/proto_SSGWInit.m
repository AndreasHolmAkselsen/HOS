clear
global M x k_cut nonLinRamp  surfaceMethod timeReached t_end H
timeReached = 0;

%% input
surfaceMethod = 'decayingConformal'; % Chalikov method
% surfaceMethod = 'Taylor';  % normal HOS

% Resolution
nx = 2^9;
M = 100; % solution order for Taylor method
relTolODE = 1e-8;

% Plot & export options
DO_EXPORT = 0;
EXPORT_MAT = 0;
exportPrefix = '';
exportPath = './figures/';
i_detailedPlot = 5; %plot contour plots of frame i. Leave empty to skip

% Wave specification
NWaves = 1;
lambda = 10;
ka = .20; % linear wave steepness

h = .1*lambda; % water depth. NB: there will be a slight differnece between H and the actual water depth when using the Chalikov method.
                % H is the water depth in the rectangular zeta-plane in this case.


% some computations...
L = NWaves*lambda;
g = 9.81;
k0 = 2*pi/lambda;
% omega = (1+.5*ka^2)*sqrt(g*k0*tanh(k0*H));
omega = sqrt(g*k0*tanh(k0*h));
T = 2*pi/omega;
c_p = 2*pi/T/k0;


% Simulation/plotting time
NT_dt = .1;
dt = NT_dt*T;
t_end = 9*dt;

    
% stability
Tramp = 0*1*T;
nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
k_cutTaylor = (M+5)*k0;
k_cut_conformal =  (M+5)*k0;%(nx*pi/L)/2;

T_init = 2*Tramp; % For the Chalikov method; the simulation will first run for time T_init with normal Taylor-HOS to generate an good initial condition.





%% Simulation

dx = L/nx;
x = (0:nx-1)'*dx - L/2;
t0 = 0;
initialStepODE = 1e-3*T;
xk0 = k0.*x;
phaseAng = 0*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);




%% init with SSGW:
N = 500;
[z,dwdz,PP] = SSGW(k0*h,ka,N);

out.hk = PP(1)*PP(2);
out.c_e = PP(4)*sqrt(g*h); % phase velocity observed from where the meam velocity at the bed is zero
out.c_s = PP(5)*sqrt(g*h); % mean flow velocity (phase velocity in frame without mean flow)
out.k = PP(2)/h;
z = [ z(N+1:end)-2*pi/(k0*h) ; z(1:N) ];
dwdz = [ dwdz(N+1:end); dwdz(1:N) ];
z = z*h;
% z = (z-z(1))*h;
dwdz = dwdz*sqrt(g*h);

n = 2*N;
z_m = .5*(z(1:n-1)+z(2:n));
dwdz0_m = .5*(dwdz(1:n-1)+dwdz(2:n))+out.c_e;
w = [0;cumsum( dwdz0_m.*diff(z))];
w = w-mean(w);

eta = interp1(real(z),imag(z),x,'linear',nan);
phiS = interp1(real(z),real(w),x,'linear',nan);


phiS0 = ka/k0.*g/omega*sin(xk0-phaseAng);

% phi = ka/k0.*g/omega*sin(xk0-phaseAng)*cosh(k*(h+z))/cosh(k*h);
u0 = ka/k0.*g/omega*k0*cos(xk0-phaseAng);
v0 =  ka/k0.*g/omega*k0*sin(xk0-phaseAng)*tanh(k0*h);
eta0 = ka/k0*(cos(xk0-phaseAng));
k_cut = k_cutTaylor;

% figure('color','w')
% subplot(311), plot(x,eta,'-',x,eta0,'--');ylabel('\eta'); grid on
% subplot(312), plot(x,phiS,'-',x,phiS0,'--');ylabel('\phi^S'); grid on
% subplot(313), plot(x,u0,'-r',x,v0,'-b',real(z),real(dwdz)+out.c_e,'--r',real(z),-imag(dwdz),'--b');ylabel('\phi^S'); grid on

if strcmp(surfaceMethod,'decayingConformal')
    
    [eta_adj,H] = initializeInitCond(x,eta,h,10);
    k_cut = k_cut_conformal;
    f = fConformal(x,eta_adj,H,inf);
    phiS_adj = interp1([x-L;x;x+L],[phiS;phiS;phiS],real(f));
    
%     figure('color','w');
%     f0 = fConformal(x,eta,H,inf);
%     fH = fConformal(x-1i*H,eta_adj,H,inf);
%     -imag(fH(1,:))
%     subplot(2,1,1); plot(x,eta,'-',real(f),imag(f),'--',real(f0),imag(f0),':','linewidth',1.5);ylabel('\eta');title('IC verification')
%     subplot(2,1,2); plot(x,phiS,'-',real(f),phiS_adj,'--','linewidth',1.5);ylabel('\phi^S');
else
    [phiS_adj,eta_adj] = deal(phiS,eta);
    [tInit,yInit] = deal([]);
    H = h;
end

% phiS_ip = phiS;
% eta_ip = eta;


tic
[t,y] = ode45(@HOSODE45 ,[t0,t_end],[phiS_adj;eta_adj],ODEoptions);
fprintf('CPU time: %gs\n',toc);

iNaN = find(isnan(y(:,1)),1,'first');
if ~isempty(iNaN), t(iNaN:end)=[]; y(iNaN:end,:)=[]; end
phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);
% interpolate to perscribed times


% t_ip = (0:dt:t_end)';
t_ip = linspace(0,t(end),10).';
nPannel = length(t_ip);
phiS_ip = interp1(t,phiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

if strcmp(surfaceMethod,'decayingConformal')
    f = fConformal(x,eta_ip,H,k_cut);
        
%     fH = fConformal(x-1i*H,eta_ip,H,k_cut);
%     -imag(fH(1,:))
%     figure, plot(x/lambda,(real(fH)-x)/lambda)
    
    for it = i_detailedPlot
        hfCont = figure('color','w'); hold on
        y0 = max(-x(end)/NWaves/2,-H);
        if y0>-H,nan_ = nan; else nan_=1; end
        
        [xi2,sig2]= ndgrid(x,linspace(y0,0,100));
        f2 = fConformal(xi2+1i*sig2,eta_ip(:,it),H,k_cut);
        haCont(1) = subplot(2,1,1);
        contour(real(f2),imag(f2),xi2,30,'r'); hold on
        contour(real(f2),imag(f2),sig2,'b');
        plot([f(:,it);nan;nan_*f2([1,end],1)],'k','linewidth',2)

        haCont(2) = subplot(2,1,2);
        kx = getKx(x);
        omega = ifft(  fft(phiS_ip(:,it)).*exp(-sig2.*kx).*2./(exp(2*kx.*H)+1).*(abs(kx)<k_cut));
        contourf(real(f2),imag(f2),real(omega),20); hold on
        contour(real(f2),imag(f2),imag(omega),20,'k')
        plot([f(:,it);nan;nan_*f2([1,end],1)],'k','linewidth',2)
        axis(haCont,'equal','off')
        
        if isfinite(H), title(haCont(1), sprintf('H[\\zeta] = %.4g, H[z] = %.4g',H,-imag(f2(1,1)))); end
        
    end

    x_ip = real(f);
%     eta_ip = imag(f); %per definition
else
    x_ip = x+0*eta_ip;
end


[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('%s Tramp%g ka=%.3g,M=%d,H=%.1f',surfaceMethod,Tramp,ka,M,H)),[],[0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),x_ip(:,i),eta_ip(:,i),'k');
    ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),nonLinRamp(t_ip(i))))
    grid(ha(i),'on');
%     axis(ha(i),'equal')
end
% linkaxes(ha)
% ylim(max(res.eta(:))*[-1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)
xlim(ha,x([1,end]));%set(ha,'XLim',x([1,end]));
set(hp(t_ip<T_init),'LineStyle','--');

    
fileName = sprintf('%s%ska%.2g_M%d_h%.2f_Nw%d_dt%.3gT_nx%d',exportPrefix,surfaceMethod,ka,M,h,NWaves,NT_dt,nx); fileName(fileName=='.')='p';
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