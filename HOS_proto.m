clear
global N dx k_cut

nx = 2^10;
g = 9.81;
N = 5; % solution order
TLin = 2.53;
NWaves = 10;
eps = .0025;


% method = 'RK4';
% method ='Euler';
method ='ODE45';
DRAW_STREAMLINE = false;

k0 = (2*pi/TLin)^2/g;
lam0 = 2*pi/k0;
L = NWaves*lam0;
dx = L/nx;
k_cut = (N+5)*k0;


x = (0:nx-1)'*dx;

% L = nx*dx;
% lam0 = L/10;
% k0 = 2*pi/lam0;
% k_cut = 5*k0;

H0 = 2*eps/k0;
% z = linspace(-3*H0,H0,50);
z = linspace(-.5*lam0,H0,50);
Psi = k0.*x;

% phiS = .5*H0.*sqrt(g/k0).*sin(k0.*x); % 3rd order Stokes for phi, not phiS!
% eps = .5*k0.*H0;
% eta = .5*H0*((1-eps^2/16)*cos(k0*x)+.5*eps*cos(2*k0*x)+3/8*eps^2*cos(3*k0*x));


%eq. 6 and 7 in HOS-memo
omega = (1+.5*eps^2).*sqrt(g*k0);
phiS = eps.*omega/k0^2*(sin(Psi)+.5*eps*sin(2*Psi) + eps^2/8*(3*sin(3*Psi)-9*sin(Psi))); 
eta = eps/k0*(cos(Psi)+.5*eps*cos(2*Psi)+3/8*eps^2*(cos(3*Psi)-cos(Psi)));
% phi0 = eps.*omega/k0^2.*(sin(Psi)+.5*eps*sin(2*Psi) + eps^2/8*(3*sin(3*Psi)-9*sin(Psi))); %.*exp(k0*z)
% [~,psi] = getStreamFunction(dx,z,fft(phi0));

%lin sol
% omega = sqrt(g*k0);
% phiS = eps.*omega/k0^2*(sin(Psi));
% eta = eps/k0*(cos(Psi));

% wave packet
% omega = sqrt(g*k0);
% % phiS = eps.*omega/k0^2*(sin(Psi)).*exp(-50*min(x/L,1-x/L).^2); 
% % eta = eps/k0*(cos(Psi)).*exp(-50*min(x/L,1-x/L).^2);
% phiS = eps.*omega/k0^2*(sin(Psi)).*exp(-10^2*(x/L-.5).^2); 
% eta = eps/k0*(cos(Psi)).*exp(-10^2*(x/L-.5).^2);


%% simulation

if strcmp(method,'ODE45')
    nPannel = 10;
    t_end = nPannel*2*pi/omega;
    [t,y] = ode45(@HOSODE45 ,[0,t_end],[phiS;eta]);
    phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);
    t_ip = linspace(0,t_end,nPannel)';
    psiS_ip = interp1(t,phiS,t_ip);
    eta_ip  = interp1(t,eta ,t_ip);
    
    
    figure('color','w','Position',[-1587 511 560 1000]); hold on;
    for iP = 1:nPannel
       subplot(nPannel,1,nPannel-iP+1), plot(x,eta_ip(iP,:),'k') ;
       ylabel(sprintf('t = %.2fs',t_ip(iP)))
    end
    return
end

dt = .001/omega;
nt = 20000;
t = 0;

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
for it = 1:nt
    switch method
        case 'Euler'
            [phiS_t,eta_t] = HOSODEeq(phiS,eta);
            phiS = phiS + phiS_t*dt;
            eta = eta + eta_t*dt;
            t = t+dt;
        case 'RK4'
            [phiS_t1,eta_t1] = HOSODEeq(phiS,eta);
            [phiS_t2,eta_t2] = HOSODEeq(phiS+.5*dt*phiS_t1, eta+.5*dt*eta_t1);
            [phiS_t3,eta_t3] = HOSODEeq(phiS+.5*dt*phiS_t2, eta+.5*dt*eta_t2);
            [phiS_t4,eta_t4] = HOSODEeq(phiS+dt*phiS_t3, eta+dt*eta_t3);
            
            phiS = phiS + (phiS_t1+2*phiS_t2+2*phiS_t3+phiS_t4)*dt/6;
            eta = eta   + (eta_t1+2*eta_t2+2*eta_t3+eta_t4)*dt/6;
            t = t+dt;
        otherwise
            error('Time integration method ''%s'' not recognised.',method)
    end
    hp.YData = eta; 
    
    if DRAW_STREAMLINE
        [~,~,~,hphi] = phiComponentsHOS(phiS,eta);
        [~,psi] = getStreamFunction(dx,z,hphi);
        hC.ZData = psi.';
    end
    drawnow
end
return



%% streamline test

[~,~,~,hphi] = phiComponentsHOS(phi0,eta);

hf = figure('color','w'); hold on;
hp = plot(x,eta,'k');
[phi,psi] = getStreamFunction(dx,z,fft(hphi));


contour(x.',z.',phi.','r');
contour(x.',z.',real(psi).','b');
axis equal
return

[W,~,~,phi,k] = phiComponentsHOS(phiS,eta);
% test that hphi(eta)=hphiS
phiS0 = 0*phiS;
for j = 1:nx
     temp = ifft(fft(phi).*exp(abs(k).*eta(j)));
     phiS0(j) = real(temp(j)); 
end
figure, plot((0:nx-1)*dx,phiS,(0:nx-1)*dx,phiS0,'--')
return

% 
% % countour plot
% z = 0:-.01:-1;
% phi = ifft(hphi.*exp(abs(k).*z),[],1);
% psi = ifft(1i*sign(k).*hphi.*exp(abs(k).*z),[],1);
% testReal( phi ) 
% testReal( psi ) 
% 
% figure('color','w'); hold on;
% contour(x.',z.',real(psi).');colorbar
% hold on
% plot(x,eta,'k')
