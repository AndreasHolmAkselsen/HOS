
addpath ./steadyState
global M x k_cut


load('./doc/figures/vortexka0p1_M5_Nw20_dt7p5T.mat')

phiS_ip = interp1(t,phiS,t_ip);
%% input


hf_S = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('S ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
hf_W = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('W ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
for iP = 2:3:nPannel
    phiS = phiS_ip(iP,:).';
    eta = eta_ip(iP,:).';
    [W_lin,W_nl,~,~,hphi,kx] = phiComponentsHOS(phiS,eta);
    W = W_lin+W_nl;
    k = abs(kx);

    % evaluate phi_x(z=eta) with manual ifft
    phi_eta = setReal(sum( hphi.*exp(k.*eta.').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
    phi_z_eta = setReal(sum( k.*hphi.*exp(k.*eta.').*exp(1i*kx.*x'), 1).'/nx,'phi_x_eta');
    phi_z_0 = setReal(ifft(k.*hphi),'phi_x_eta');
    phi_0 = setReal(ifft(hphi),'phi_x_eta');
    
    
    figure(hf_S)
    subplot(nPannel,1,nPannel-iP+1);
    plot(x,W,'k',x,phi_z_eta,'--r',x,phi_z_0,':b','linewidth',1);    legend('W','\phi_z(z=\eta)','\phi_z(z=0)');
        set(gca,'XTick',[]);
    box off; grid on;
    
    figure(hf_W)
    subplot(nPannel,1,nPannel-iP+1);
    plot(x,phiS.','k',x,phi_eta,'--r',x,phi_0,':b','linewidth',1);    legend('\phi^S','\phi(z=\eta)','\phi(z=0)');
    set(gca,'XTick',[]);
    box off; grid on;
end

