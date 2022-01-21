clear
addpath ./steadyState
global M x k_cut
g = 9.81;

%% input
nx = 2^10;
M = 10; % solution order

L = 100;
dk = 2*pi/L; kMax = dk*floor((nx-1)/2);
k_cut = kMax/2;

dx = L/nx;
x = (0:nx-1)'*dx;

eta = L*.005*exp(-(x-L/2).^2/(.1*L)^2);
phiS = L*.01*exp(-(x-.4*L).^2/(.12*L)^2).*cos(10*2*pi*x/L);

% [phiS_x,eta_x,phiS,hphi,kx] = phiComponentsHOS_steadyState(W,eta);

[W_lin,W_nl,phiS_x,eta_x,hphi,kx] = phiComponentsHOS(phiS,eta);
W = W_lin+W_nl;

% evaluate phi_x(z=eta) with manual ifft
k = abs(kx);
phi_eta = setReal(sum( hphi.*exp(k.*eta').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
phi_z_eta = setReal(sum( k.*hphi.*exp(k.*eta').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
phi_z_0 = setReal(ifft(k.*hphi),'phi_x_eta');
phi_0 = setReal(ifft(hphi),'phi_x_eta');

% phi_x_eta = setReal(sum( 1i*kx.*hphi.*exp(k.*eta').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
% phi_x_0 = setReal(ifft(1i*kx.*hphi),'phi_x_eta');

hf = figure('color','w','position',[-1919   401  1280   603]);
subplot(3,1,1);plot(x,eta,'k','linewidth',1);    ylabel('\eta');
subplot(3,1,2);plot(x,W,'k',x,phi_z_eta,'--r',x,phi_z_0,':b','linewidth',1);    legend('W','\phi_z(z=\eta)','\phi_z(z=0)');
subplot(3,1,3);plot(x,phiS,'k',x,phi_eta,'--r',x,phi_0,':b','linewidth',1);    legend('\phi^S','\phi(z=\eta)','\phi(z=0)');