
addpath ./steadyState
global M x k_cut

% matFilePath = './doc/figures/vortexka0p2_M5_Nw20_dt1T';
matFilePath = './doc/figures/vortexka0p1_M5_Nw20_dt7p5T';
load([matFilePath,'.mat']);
phiS_ip = interp1(t,phiS,t_ip);

%% input

for iP = 8:10
    hf = figure('color','w','Position',[649   359   522   446],'name',sprintf('ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
    phiS = phiS_ip(iP,:).';
    eta = eta_ip(iP,:).';
    t = t_ip(iP);
    [W_lin,W_nl,~,~,hphi,kx] = phiComponentsHOS(phiS,eta);
    W = W_lin+W_nl;
    k = abs(kx);

    % evaluate phi_x(z=eta) with manual ifft
    phi_eta = setReal(sum( hphi.*exp(k.*eta.').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
    phi_z_eta = setReal(sum( k.*hphi.*exp(k.*eta.').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
    phi_z_0 = setReal(ifft(k.*hphi),'phi_x_eta');
    phi_0 = setReal(ifft(hphi),'phi_x_eta');


    ha(1)=subplot(3,1,1);plot(x,eta,'k','linewidth',1);    legend('\eta','fontsize',10);  grid on; title(sprintf('Time: %.1fs',t))
    ha(2)=subplot(3,1,2);plot(x,phiS,'k',x,phi_eta,'--r',x,phi_0,':b','linewidth',1);    legend({'\phi^S','\phi(z=\eta)','\phi(z=0)'},'fontsize',10); grid on;
    ha(3)=subplot(3,1,3);plot(x,W,'k',x,phi_z_eta,'--r',x,phi_z_0,':b','linewidth',1);    legend({'W','\phi_z(z=\eta)','\phi_z(z=0)'},'fontsize',10); grid on;

    set(ha(1:2),'XTick',[]);
    axis(ha,'tight');
    set(ha,'Box','off','XLim',[.45,.55]*L);
    linkaxes(ha,'x')

    if ~isfolder([matFilePath,'_Taylor']),mkdir([matFilePath,'_Taylor']); end
    export_fig([matFilePath,'_Taylor/t_',num2str(t,'%.0f')],'-pdf')
end