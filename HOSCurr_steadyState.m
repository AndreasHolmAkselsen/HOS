clear
addpath ./steadyState
global M dx k_cut 
g = 9.81;

%% input
nx = 2^10;
PLOT_CURRENT = false; 
PLOT_FINAL_VALOCITY_FIELD = true;

L = 100;
dk = 2*pi/L; kMax = dk*floor((nx-1)/2);


M = 10; % solution order
nIt = 20;
relax = 1;
k_cut = kMax/2;

nMirror = 3;

% % similar to basin
% ID_j   = [       1,     -1  ];% 1: Source/sink, -1: vortex  
% zeta_j = [-.05-.3i,  .1-.1i ]*L;% object centre
% F_j    = [       3,     1   ];% object strength


% single vortex
ID_j   = [  -1, -1 ];% 1: Source/sink, -1: vortex  
zeta_j = [.3-.05i, .7-.05i ]*L;% object centre
Fr_j   = .1*[ 1, .7];% object strength, in Froude units

dx = L/nx;
x = (0:nx-1)'*dx;

Fr_j = shiftdim(Fr_j,-1); zeta_j = shiftdim(zeta_j,-1); ID_j = shiftdim(ID_j,-1);
zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);
A_j = Fr_j.*sqrt(g).*abs(imag(zeta_j)).^1.5;
f = @(zeta) sum(A_j.*sqrt(ID_j).*(ID_j.*log(zeta-zeta_j) + log(zeta-conj(zeta_j))),3:4);
df = @(zeta) sum(A_j.*sqrt(ID_j).*(ID_j./(zeta-zeta_j) + 1./(zeta-conj(zeta_j))),3:4);

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



%% simulation

[eta,phiS,W,hphi,kx] =  HOSODEeqCurr_steadyState(x,df,nIt,relax);

if PLOT_FINAL_VALOCITY_FIELD
    z = linspace(-2*max(abs(eta)),max(eta),500);
%     z = linspace(.9*max(imag(zeta_j(:,:,1))),max(eta),500);
    
    [phi,psi] = getStreamFunction(dx,z,hphi);
    psiFull = psi+imag(f(x+1i*z));
    psiFull(z>eta) = nan;
    
    figure('color','w');
    contour(x',z',psiFull',20)
    hold on; 
    plot(x,eta,'r','linewidth',1);    ylabel('\eta');
    
    % evaluate phi_x(z=eta) with manual ifft
    k = abs(kx);
    phi_x_eta = setReal(sum( 1i*kx.*hphi.*exp(k.*eta').*exp(1i*kx.*x').*(k<k_cut), 1).'/nx,'phi_x_eta');
    phi_z_eta = setReal(sum( k.*hphi.*exp(k.*eta').*exp(1i*kx.*x').*(k<k_cut), 1).'/nx,'phi_x_eta');
    phi_x_0 = setReal(ifft(1i*kx.*hphi.*(k<k_cut)),'phi_x_eta');
    phi_z_0 = setReal(ifft(k.*hphi.*(k<k_cut)),'phi_x_eta');
    eta_x = setReal(ifft(1i*kx.*fft(eta).*(k<k_cut)),'eta');
    
    U_eta = conj(df(x+1i*eta));
    u_eta = real(U_eta)+phi_x_eta;
    w_eta =  imag(U_eta)+phi_z_eta;
    wBerboulli = .5*max(abs(df(x)))^2;
    bernoulli = .5*u_eta.^2 + .5*w_eta.^2 + g*eta;
    
    kinBC = u_eta.*eta_x - w_eta;
    wKin = max(abs(df(x)));
    
    figure('color','w','position',[-1919   401  1280   603]);    
    subplot(5,1,1);plot(x,eta,'k','linewidth',1);    ylabel('\eta');
    title("M = "+M+", nIt = "+nIt+", relax = "+relax)
    subplot(5,1,2);plot(x,bernoulli/wBerboulli,'k','linewidth',1);   
    ylabel('$\frac{\mathrm{Bernoulli}}{\frac12|U|_\mathrm{max}^2}$','interpreter','latex','fontsize',14);
    subplot(5,1,3);plot(x,kinBC/wKin,'k','linewidth',1);    
    ylabel('$\frac{\mathrm{kin.\, BC}}{|U|_\mathrm{max}}$','interpreter','latex','fontsize',14);
    subplot(5,1,4);plot(x,W,'k',x,phi_z_eta,'--r',x,phi_z_0,':b','linewidth',1);    legend('W','\phi_z(z=\eta)','\phi_z(z=0)');
    subplot(5,1,5);plot(x,phi_x_eta,'--r',x,phi_x_0,':b','linewidth',1);    legend('\phi_x(z=\eta)','\phi_x(z=0)');
end


