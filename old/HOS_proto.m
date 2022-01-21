clear
global M dx k_cut Tramp nRamp HOSODEsyst
g = 9.81;

%% input
nx = 2^10;
M = 5; % solution order

NWaves = 10;
lambda = 10;
ka = .28;

% some computations...
k0 = 2*pi/lambda;
T = 2*pi/((1+.5*ka^2)*sqrt(g*k0));

L = NWaves*lambda;

dt = T;
t_end = 9*dt;


% Ramp: wNl = 1-exp(-(t/Tramp)^nRamp);
Tramp = 2*T;
nRamp = 2;
k_cut = (M+5)*k0;

method ='ODE45'; % {'ODE45','RK4','Euler'}
IC = 'linearWave'; % {'linearWave','Stokes3','wavePacket'}

% if IC='wavePacket':
packetWidth = .1*L;
x0 = 2/5*L;

relTolODE = 1e-8;
DRAW_STREAMLINE = false; %outdated?

%% code

dx = L/nx;
x = (0:nx-1)'*dx;

H0 = 2*ka/k0;
if DRAW_STREAMLINE, z = linspace(-.5*lam0,H0,50); end
xk0 = k0.*x;

omegaLin = sqrt(g*k0);
switch IC
    case 'Stokes3'
%         eq. 6 and 7 in HOS-memo
        omega = (1+.5*ka^2).*sqrt(g*k0);
        phiS = ka.*omega/k0^2*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0)));
        eta = ka/k0*(cos(xk0)+.5*ka*cos(2*xk0)+3/8*ka^2*(cos(3*xk0)-cos(xk0)));
        % phi0 = ka.*omega/k0^2.*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0))); %.*exp(k0*z)
        % [~,psi] = getStreamFunction(dx,z,fft(phi0));
    case 'linearWave'
        phiS = ka.*omegaLin/k0^2*(sin(xk0));
        eta = ka/k0*(cos(xk0));
    case 'wavePacket'
        packet = exp(-(min(abs(x-x0),L-abs(x-x0))/packetWidth).^2);
        phiS = ka.*omegaLin/k0^2*(sin(xk0)).*packet;
        eta = ka/k0*(cos(xk0)).*packet; 
    otherwise
        
end

%% simulation

if strcmp(method,'ODE45')
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
    for iP = 1:nPannel
       subplot(nPannel,1,nPannel-iP+1), plot(x,eta_ip(iP,:),'k');hold on
       ylabel(sprintf('t = %.2fs',t_ip(iP)))
    end
    return
end


hf = figure('color','w','Position',[-1587 511 560 400]); hold on;
hp = plot(x,eta,'k');
if DRAW_STREAMLINE
    [~,hC] = contour(x.',z.',0*(x.*z).');
    ylim([z(1),1.2*max(eta)]);
else
    ylim(1.2*[min(eta),max(eta)]);
end
% axis equal manual
drawnow;


t = 0;
while t <= t_end
    switch method
        case 'Euler'
            [phiS_t,eta_t] = HOSODEeq(t,phiS,eta);
            phiS = phiS + phiS_t*dt;
            eta = eta + eta_t*dt;
        case 'RK4'
            [phiS_t1,eta_t1] = HOSODEeq(t,phiS,eta);
            [phiS_t2,eta_t2] = HOSODEeq(t,phiS+.5*dt*phiS_t1, eta+.5*dt*eta_t1);
            [phiS_t3,eta_t3] = HOSODEeq(t,phiS+.5*dt*phiS_t2, eta+.5*dt*eta_t2);
            [phiS_t4,eta_t4] = HOSODEeq(t,phiS+dt*phiS_t3, eta+dt*eta_t3);
            
            phiS = phiS + (phiS_t1+2*phiS_t2+2*phiS_t3+phiS_t4)*dt/6;
            eta = eta   + (eta_t1+2*eta_t2+2*eta_t3+eta_t4)*dt/6;
        otherwise
            error('Time integration method ''%s'' not recognised.',method)
    end
    hp.YData = eta; 
    t = t+dt;
    if DRAW_STREAMLINE
        [~,~,~,~,hphi] = phiComponentsHOS(phiS,eta);
        [~,psi] = getStreamFunction(dx,z,hphi);
        hC.ZData = psi.';
    end
    drawnow
end
return



%% streamline test

[~,~,~,~, hphi] = phiComponentsHOS(phi0,eta);

hf = figure('color','w'); hold on;
hp = plot(x,eta,'k');
[phi,psi] = getStreamFunction(dx,z,fft(hphi));


contour(x.',z.',phi.','r');
contour(x.',z.',real(psi).','b');
axis equal
return

[~,~,~,~,phi,k] = phiComponentsHOS(phiS,eta);
% test that hphi(eta)=hphiS
phiS0 = 0*phiS;
for j = 1:nx
     temp = ifft(fft(phi).*exp(abs(k).*eta(j)));
     phiS0(j) = real(temp(j)); 
end
figure, plot((0:nx-1)*dx,phiS,(0:nx-1)*dx,phiS0,'--')
return

