clear
global M x k_cut nonLinRamp 

%% input
% Resolution
nx = 2^9;
M = 10; % solution order
relTolODE = 1e-8;

% Plot & export options
DO_EXPORT = false;
exportPrefix = 'AHA_';
exportPath = './HOS_SFo_curr/figures/';

% Wave specification
NWaves = 1;
lambda = 1;
ka = .35;

L = NWaves*lambda;
initialCondition = 'linearWave'; % {'linearWave','Stokes3','wavePacket'} 
% if initialCondition='wavePacket':
% packetWidth = .1*L;
% x0 = 2/5*L;

% some computations...
g = 9.81;
k0 = 2*pi/lambda;
omega = (1+.5*ka^2)*sqrt(g*k0);
T = 2*pi/omega;
c_p = 2*pi/T/k0;

% Simulation/plotting time
NT_dt = 1;
dt = NT_dt*T;
t_end = 9*dt;

% stability
Tramp = 1*T;
nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
k_cut = (M+5)*k0;

initialStepODE = 1e-3*T;


%% Plot background current
dx = L/nx;
x = (0:nx-1)'*dx;

%% Simulation
xk0 = k0.*x;
phaseAng = 30*pi/180;
switch initialCondition
    case 'Stokes3'
%         eq. 6 and 7 in HOS-memo
        phiS = ka.*omega/k0^2*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0)));
        eta = ka/k0*(cos(xk0)+.5*ka*cos(2*xk0)+3/8*ka^2*(cos(3*xk0)-cos(xk0)));
        % phi0 = ka.*omega/k0^2.*(sin(xk0)+.5*ka*sin(2*xk0) + ka^2/8*(3*sin(3*xk0)-9*sin(xk0))); %.*exp(k0*z)
        % [~,psi] = getStreamFunction(dx,z,fft(phi0));
    case 'linearWave'
        phiS = ka.*sqrt(g*k0)/k0^2*(sin(xk0-phaseAng));
        eta = ka/k0*(cos(xk0-phaseAng));
%         phiS = ka.*sqrt(g*k0)/k0^2*(sin(xk0-phaseAng)+0*sin(2*xk0-100*rand)+0*sin(3*xk0-100*rand)+.0*sin(4*xk0-100*rand));
%         eta = ka/k0*(               0*cos(xk0-phaseAng)+0*cos(2*xk0-100*rand)+.0*cos(3*xk0-100*rand)+.1*cos(4*xk0-100*rand));
    case 'wavePacket'
        packet = exp(-(min(abs(x-x0),L-abs(x-x0))/packetWidth).^2);
        phiS = ka.*sqrt(g*k0)/k0^2*(sin(xk0)).*packet;
        eta = ka/k0*(cos(xk0)).*packet; 
    otherwise
        
end


[W_lin,W_nl,phiS_x,eta_x,hphi,kx] = phiComponentsHOS(phiS,eta,inf,M);
        
[U,eta_x,homega,kxp] = phiComponentsConformal(phiS,eta);
psiSConf = real(sum(homega.*exp(-1i*kxp.*x.'),1))';
WHOS = W_lin+W_nl;
WConf = imag(U);

% compute dw/dz using finite difference and compare to dw/dz = omega'/f'.
k = abs(kx);
heta = fft(eta)/nx.*(k<k_cut);
eps = .05*dx*exp(1i*pi/3); % some random phase.
xi = x;
zz1 = xi - eps;
zz2 = xi + eps;
z1 = zz1 + 1i*sum(heta.'.*exp(1i*kx.'.*zz1),2);
z2 = zz2 + 1i*sum(heta.'.*exp(1i*kx.'.*zz2),2);
dw = ( sum(homega.'.*exp(-1i*kxp.'.*zz2),2)-sum(homega.'.*exp(-1i*kxp.'.*zz1),2) )...
     ./( z2 - z1  );

figure,
subplot(2,1,1); plot(x,phiS,x,psiSConf,'--','linewidth',1.5);ylabel('\phi^S');title('conformal test');
subplot(2,1,2); plot(x,-imag(dw),x,WConf,'--','linewidth',1.5);ylabel('W');


% compute phi and w at z=x+i eta from the Taylor potential field hphi. Compare to phiS and W.
phiS0 = real(sum(hphi.*exp(1i*kx.*x.'+abs(kx).*eta.'),1).'/nx);
W0 = real(sum(abs(kx).*hphi.*exp(1i*kx.*x.'+abs(kx).*eta.'),1).'/nx);
U0 = real(sum(1i*kx.*hphi.*exp(1i*kx.*x.'+abs(kx).*eta.'),1).'/nx);

% phiyy0 = real(sum(abs(kx).^2.*hphi.*exp(1i*kx.*x.'+abs(kx).*eta.'),1).'/nx);
% phixx0 = real(sum(-kx.^2.*hphi.*exp(1i*kx.*x.'+abs(kx).*eta.'),1).'/nx);
% figure, plot(x,phiyy0+phixx0)

figure,
subplot(2,1,1);plot(x,phiS,x,phiS0,'--','linewidth',1.5);ylabel('\phi^S');title('Taylor test');
subplot(2,1,2);plot(x,WHOS,x,W0,'--','linewidth',1.5);ylabel('W');



figure, 
subplot(2,1,1);plot(x,WHOS,'',x,WConf,'--',x,W_lin,':','linewidth',1.5);title('Taylor vs conformal, main')
subplot(2,1,2);plot(x,W0,'',x,-imag(dw),'--','linewidth',1.5);title('Taylor vs conformal, ref.')

figure;
subplot(2,1,1);plot(x,U0,'-',x, real(U),'--',x,real(dw),':','linewidth',1.5);ylabel('u');
subplot(2,1,2);plot(x,W0,'-',x, imag(U),'--',x,-imag(dw),':','linewidth',1.5);ylabel('w');
% subplot(2,1,3);plot(x,sqrt(phix0.^2 + W0.^2),'-',x, abs(U),'--',x, abs(dw),':','linewidth',1.5);ylabel('Abs U');

% [y,sigma] = deal(linspace(2*min(eta),max(eta),100));
[y,sigma] = deal(2*min(eta):dx:max(eta));
zeta = xi+1i*sigma;
omega = sum(shiftdim(homega,-2).*exp(-1i*shiftdim(kxp,-2).*zeta),3);
z = zeta + 1i.*sum(shiftdim(heta,-2).*exp(1i*shiftdim(kx,-2).*zeta),3);
eta_xi = interp1(x,eta,real(z),'linear','extrap');
omega(imag(z)>eta_xi+.001) = nan;

% phi = real(sum(shiftdim(hphi,-2).*exp(1i*shiftdim(kx,-2).*x + shiftdim(k,-2).*z ),3))/nx;
phi = real(ifft(hphi.*exp(k.*y),[],1));
phi( y > eta ) = nan;

% figure
% subplot(1,2,1);contourf(real(z),imag(z),real(omega),30,'LineStyle','none');hold on; 
% surface(x.*[1,1],eta.*[1,1],0*[x,x],[phiS,phiS],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2); colorbar
%     
% subplot(1,2,2);contourf(x+0*y,y+0*x,phi,30,'LineStyle','none');hold on; 
% surface(x.*[1,1],eta.*[1,1],0*[x,x],[phiS,phiS],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2); colorbar    
    
figure
xc = .25*max(x);
omegaNan = omega;
omegaNan(real(z)>xc) = nan;
cl = linspace(min(phiS),max(phiS),30);
[~,hc1]=contourf(real(z),imag(z),real(omegaNan),cl,'LineStyle','none');hold on; 
[~,hc2]=contourf(x(x>xc)+0*y,y+0*x(x>xc),phi(x>xc,:),cl,'LineStyle','none');
surface(x.*[1,1],eta.*[1,1],0*[x,x],[phiS,phiS],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);         
    
    
figure
yc = 1.2*min(eta);
omegaNan = omega;
omegaNan(imag(z)>yc) = nan;
cl = linspace(min(phiS),max(phiS),30);
[~,hc1]=contourf(real(z),imag(z),real(omegaNan),cl,'LineStyle','none');hold on; 
[~,hc2]=contourf(x+0*y(y>yc),y(y>yc)+0*x,phi(:,y>yc),cl,'LineStyle','none');
surface(x.*[1,1],eta.*[1,1],0*[x,x],[phiS,phiS],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);             
    
    
zzp = 4.26 + (-2:.01:0)*1i;
zp = zzp + 1i.*sum(shiftdim(heta,-2).*exp(1i*shiftdim(kx,-2).*zzp),3);
phip =  real(sum(shiftdim(hphi,-2).*exp(1i*shiftdim(kx,-2).*real(zp) + shiftdim(k,-2).*imag(zp) ),3))/nx;
omegap = sum(shiftdim(homega,-2).*exp(-1i*shiftdim(kxp,-2).*zzp),3);
figure,
subplot(2,1,1);plot(phip,imag(zp),'-',real(omegap),imag(zp),'--');ylabel('y');xlabel('\phi');


zzp = 0:.001:max(x)  - 0.23i;
zp = zzp + 1i.*sum(shiftdim(heta,-2).*exp(1i*shiftdim(kx,-2).*zzp),3);
phip =  real(sum(shiftdim(hphi,-2).*exp(1i*shiftdim(kx,-2).*real(zp) + shiftdim(k,-2).*imag(zp) ),3))/nx;
omegap = sum(shiftdim(homega,-2).*exp(-1i*shiftdim(kxp,-2).*zzp),3);
subplot(2,1,2); plot(real(zp),phip,'-',real(zp),real(omegap),'--');ylabel('\phi');xlabel('x');



% figure, hold on;
% for i=1:8
% zzp = (0:.01:max(x))'  - i*1.2i*max(eta);
% zp = zzp + 1i.*sum(shiftdim(heta,-2).*exp(1i*shiftdim(kx,-2).*zzp),3);
% omegap = sum(shiftdim(homega,-2).*exp(-1i*shiftdim(kxp,-2).*zzp),3);
%  plot(real(zp),imag(zp))
% % figure; plot(real(zp),real([omegap]))
% % surface(real(zp).*[1,1],imag(zp).*[1,1],real(zp)*[0,0],real([omegap,omegap]),...
% %         'facecol','no',...
% %         'edgecol','interp',...
% %         'linew',2);   colorbar
% end