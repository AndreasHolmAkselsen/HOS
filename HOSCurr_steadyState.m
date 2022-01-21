clear
addpath ./steadyState
global M dx k_cut 
g = 9.81;

%% input
nx = 2^10;
PLOT_CURRENT = true; 
PLOT_FINAL_VALOCITY_FIELD = true;
DO_EXPORT = true;
exportPrefix = 'steady_2vortex_';
exportPath = './doc/figures/';

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
zeta_j = [.3-.05i, .7-.05i ]*L;% object centre
Fr_j   = [ .2i, .15i];% object strength, in Froude units

dx = L/nx;
x = (0:nx-1)'*dx;

Fr_j = shiftdim(Fr_j,-1); zeta_j = shiftdim(zeta_j,-1); %ID_j = shiftdim(ID_j,-1);
zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);
% A_j = Fr_j.*sqrt(g).*abs(imag(zeta_j)).^1.5;
% f = @(zeta) sum(A_j.*sqrt(ID_j).*(ID_j.*log(zeta-zeta_j) + log(zeta-conj(zeta_j))),3:4);
% df = @(zeta) sum(A_j.*sqrt(ID_j).*(ID_j./(zeta-zeta_j) + 1./(zeta-conj(zeta_j))),3:4);

A_j = .5*Fr_j.*sqrt(g).*abs(imag(zeta_j)).^1.5;
f  = @(zeta) sum(A_j.*log(zeta-zeta_j) + conj(A_j.*log(conj(zeta)-zeta_j)),3:4);
df = @(zeta) sum(A_j./(zeta-zeta_j) + conj(A_j./(conj(zeta)-zeta_j)),3:4);

if PLOT_CURRENT 
    z = linspace(-.3*L,.3*L,100);
    hf_c = figure('color','w'); ha = gca;
    % plot velocity intensity |U|
    absU = abs(df(x+1i*z));
    ULim = 2.5*max(abs(df(x)));
    hIm = imagesc(x',z',absU',[0,ULim]); 
%     colorbar
    hIm.AlphaData = (absU<ULim)';
    % plot streamlines
    hold on, axis equal xy
    contour(x',z',imag(f(x+1i*z))',20,'k');
    plot(ha.XLim,[0,0],'k')
    ylim([z(1),0])
    
    % add quiver plot
    zeta_ip = linspace(x(1),x(end),10) + 1i*linspace(z(1),z(end),10)';
    df_ip = df(zeta_ip);
    df_ip(abs(df_ip)>ULim) = nan;
    quiver(real(zeta_ip),imag(zeta_ip),real(df_ip),-imag(df_ip),'r');
    
    ha.Visible='off';
    
%     % plot horizontal velocity at z=0
%     figure('color','w');
%     plot(x,real(df(x)))
%     title(sprintf('|\\phi_x^{(1)}| = %.2gm/s, c_p = %.2gm/s',ka*sqrt(g/k0),c_p));
    
    fprintf('max |U(0)| = %.3gm/s\n',max(abs(df(x))))
    drawnow
    
    if DO_EXPORT
        fileName = sprintf('curr_%sM%d',exportPrefix,M); fileName(fileName=='.')='p';
        savefig(hf_c,[exportPath,'/',fileName]);
        export_fig(hf_c,[exportPath,'/',fileName],'-pdf');
    end
    
end


%% simulation

[eta,phiS,W,hphi,kx] =  HOSODEeqCurr_steadyState(x,df,nIt,relax);
% [eta_lin,phi0_lin,hphi_lin] =  HOSODEeqCurr_linearSteadyState(x,df,nIt);


if PLOT_FINAL_VALOCITY_FIELD
    z = linspace(-2*max(abs(eta)),max(eta),500);
%     z = linspace(.9*max(imag(zeta_j(:,:,1))),max(eta),500);
    
    [phi,psi] = getStreamFunction(dx,z,hphi);
    psiFull = psi+imag(f(x+1i*z));
    psiFull(z>eta) = nan;
    
    hf_sl = figure('color','w','Position',[-1417 816 560 98]);
    contour(x',z',psiFull',20)
    colormap(flipud(colormap))
    hold on; 
    plot(x,eta,'k','linewidth',1);    
%     ylabel('$\eta$','interpreter','latex','fontsize',14);
    set(gca,'box','off');%'XTick',[],
    
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
    
    hf_err = figure('color','w','position',[ -1686 698 574 253]);
    ha(1)=subplot(2,1,1);plot(x,kinBC/wKin,'k','linewidth',1);
    ylabel('$\frac{\mathrm{kin.\, BC}}{|U|_\mathrm{max}}$','interpreter','latex','fontsize',14);grid on;
    ha(2)=subplot(2,1,2);plot(x,bernoulli/wBerboulli,'k','linewidth',1);   grid on;
    ylabel('$\frac{\mathrm{Bernoulli}}{\frac12|U|_\mathrm{max}^2}$','interpreter','latex','fontsize',14);
    set(ha,'XTick',[])
    
    
%         hf_err = figure('color','w','position',[-1686 592 574 359]);%[-1919   401  1280   603]);    
%     ha(1)=subplot(3,1,1);plot(x,eta,'k','linewidth',1);    ylabel('$\eta$','interpreter','latex','fontsize',14); grid on;
% %     title("M = "+M+", nIt = "+nIt+", relax = "+relax)
%     ha(2)=subplot(3,1,2);plot(x,bernoulli/wBerboulli,'k','linewidth',1);   grid on;
%     ylabel('$\frac{\mathrm{Bernoulli}}{\frac12|U|_\mathrm{max}^2}$','interpreter','latex','fontsize',14);
%     ha(3)=subplot(3,1,3);plot(x,kinBC/wKin,'k','linewidth',1);    
%     ylabel('$\frac{\mathrm{kin.\, BC}}{|U|_\mathrm{max}}$','interpreter','latex','fontsize',14);grid on;
% %     subplot(5,1,4);plot(x,W,'k',x,phi_z_eta,'--r',x,phi_z_0,':b','linewidth',1);    legend('W','\phi_z(z=\eta)','\phi_z(z=0)');
% %     subplot(5,1,5);plot(x,phi_x_eta,'--r',x,phi_x_0,':b','linewidth',1);    legend('\phi_x(z=\eta)','\phi_x(z=0)');
%     set(ha,'XTick',[])
    

    if DO_EXPORT
        fileName = sprintf('%sM%d',exportPrefix,M); fileName(fileName=='.')='p';
        savefig(hf_err,[exportPath,'/',fileName]);
        export_fig(hf_err,[exportPath,'/',fileName],'-pdf');
        
        savefig(hf_sl,[exportPath,'/sl_',fileName]);
        export_fig(hf_sl,[exportPath,'/sl_',fileName],'-pdf');
    end
    
    
%     % Error in linear solution
% %     hphi_lin = 0*hphi_lin;
%     phi_x_eta_lin = setReal(sum( 1i*kx.*hphi_lin.*exp(k.*eta_lin').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
%     phi_z_eta_lin = setReal(sum( k.*hphi_lin.*exp(k.*eta_lin').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
%     eta_x_lin = setReal(ifft(1i*kx.*fft(eta_lin)),'eta');
%     U_eta_lin = conj(df(x+1i*eta_lin));
%     u_eta_lin = real(U_eta_lin)+phi_x_eta_lin;
%     w_eta_lin =  imag(U_eta_lin)+phi_z_eta_lin;
%     bernoulli_lin = .5*u_eta_lin.^2 + .5*w_eta_lin.^2 + g*eta_lin;
%     kinBC_lin = u_eta_lin.*eta_x_lin - w_eta_lin;
% 
%     
%     hf_linErr = figure('color','w','position',[-1919   401  1280   603]);    
%     subplot(3,1,1);plot(x,eta_lin,'k','linewidth',1);    ylabel('\eta');
%     title('Linear solution')
%     subplot(3,1,2);plot(x,bernoulli_lin/wBerboulli,'k','linewidth',1);   
%     ylabel('$\frac{\mathrm{Bernoulli}}{\frac12|U|_\mathrm{max}^2}$','interpreter','latex','fontsize',14);
%     subplot(3,1,3);plot(x,kinBC_lin/wKin,'k','linewidth',1);    
%     ylabel('$\frac{\mathrm{kin.\, BC}}{|U|_\mathrm{max}}$','interpreter','latex','fontsize',14);



end


