clear
global M dx x k_cut Tramp nRamp HOSODEsyst df
g = 9.81;

%% input
nx = 2^10;
M = 5; % solution order
PLOT_CURRENT = true; 
PLOT_FINAL_VALOCITY_FIELD = true;

NWaves = 20;
lambda = 10;
ka = 0;%.2;

% some computations...
k0 = 2*pi/lambda;
T = 2*pi/((1+.5*ka^2)*sqrt(g*k0));
c_p = 2*pi/T/k0;

L = NWaves*lambda;

dt = T;
t_end = 9*dt;


% Ramp: wNl = 1-exp(-(t/Tramp)^nRamp);
Tramp = 2*T;
nRamp = 2;
k_cut = (M+5)*k0;

IC = 'linearWave'; % {'linearWave','Stokes3','wavePacket'}

% if IC='wavePacket':
packetWidth = .1*L;
x0 = 2/5*L;

relTolODE = 1e-8;
nMirror = 3;

% % similar to basin
% ID_j   = [       1,     -1  ];% 1: Source/sink, -1: vortex  
% zeta_j = [-.05-.3i,  .1-.1i ]*L;% object centre
% F_j    = c_p*20*[       3,     1   ];% object strength


% % single vortex
ID_j   = [      -1 ];% 1: Source/sink, -1: vortex  
zeta_j = [.5-.075i  ]*L;% object centre
F_j    = c_p*20*[    1  ];% object strength

% single vortex
% ID_j   = [      1 ,1];% 1: Source/sink, -1: vortex  
% zeta_j = [.15-.7i ,.5-.7i  ]*L;% object centre
% F_j    = c_p*20*[    1,-1  ];% object strength



dx = L/nx;
x = (0:nx-1)'*dx;

% x = (-nx:2*nx-1)'*dx;


F_j = shiftdim(F_j,-1); zeta_j = shiftdim(zeta_j,-1); ID_j = shiftdim(ID_j,-1);
zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);

f = @(zeta) sum(F_j.*sqrt(ID_j)./(2*pi).*(ID_j.*log(zeta-zeta_j) + log(zeta-conj(zeta_j))),3:4);
df = @(zeta) sum(F_j.*sqrt(ID_j)./(2*pi).*(ID_j./(zeta-zeta_j) + 1./(zeta-conj(zeta_j))),3:4);

if PLOT_CURRENT
    z = linspace(-.3*L,.3*L,100);
    figure('color','w'); ha = gca;
    % plot velocity intensity |U|
    absU = abs(df(x+1i*z));
    ULim = 2.5*max(abs(df(x)));
    hIm = imagesc(x',z',absU',[0,ULim]); colorbar
    hIm.AlphaData = (absU<ULim)';
    % plot streamlines
    hold on, axis equal
    contour(x',z',imag(f(x+1i*z))',20,'k');
    
    plot(ha.XLim,[0,0],'k')
    
    % add quiver plot
    zeta_ip = linspace(x(1),x(end),10) + 1i*linspace(z(1),z(end),10)';
    df_ip = df(zeta_ip);
    quiver(real(zeta_ip),imag(zeta_ip),real(df_ip),-imag(df_ip),'r');
    
    % plot horizontal velocity at z=0
    figure('color','w');
    plot(x,real(df(x)))
    title(sprintf('|\\phi_x^{(1)}| = %.2gm/s, c_p = %.2gm/s',ka*sqrt(g/k0),c_p));
    drawnow
end

%% code
xk0 = k0.*x;
switch IC
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

%% simulation
HOSODEsyst = @HOSODEeqCurr;
ODEoptions = odeset('RelTol',relTolODE);
tic
[t,y] = ode45(@HOSODE45 ,[0,t_end],[phiS;eta],ODEoptions);
fprintf('CPU time (AHA): %gs\n',toc);
phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);
t_ip = (0:dt:t_end)';
nPannel = length(t_ip);
psiS_ip = interp1(t,phiS,t_ip);
eta_ip  = interp1(t,eta ,t_ip);


figure('color','w','Position',[-1587 511 560 1000]);
wNc=max(0,1-exp(-(t_ip/Tramp).^nRamp));
for iP = 1:nPannel
    subplot(nPannel,1,nPannel-iP+1), plot(x,eta_ip(iP,:),'k');hold on
    ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),wNc(iP)))
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
