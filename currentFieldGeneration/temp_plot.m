cases = {
    'T1p0_ka0p2_M5_Nw69_dt5T_L130'
    'T1p5_ka0p2_M5_Nw32_dt5T_L130'
    'T2p0_ka0p2_M5_Nw19_dt5T_L130'
    'T2p5_ka0p2_M5_Nw12_dt5T_L130'
    'T2p9_ka0p2_M5_Nw9_dt5T_L130'
    'T3p6_ka0p1_M5_Nw6_dt5T_L130'
    'T4p0_ka0p073_M5_Nw5_dt5T_L130'
    'T4p5_ka0p049_M5_Nw4_dt5T_L130'
    };

filePath = '.\doc\figures\basin\';
DO_EXPORT = false;
PLOT_WITH_SUBPLOT = true;

for i = 1:length(cases)
    load([filePath,cases{i},'.mat']);
    
    iEnd = find(any(isnan(eta),2),1,'first');
    if ~isempty(iEnd)
        t_end = t(iEnd);
    else
        t_end = t(end);
    end
    dt = t_end/9;
    % interpolate to perscribed times
    t_ip = (0:dt:t_end)';
    nPannel = length(t_ip);
    % phiS_ip = interp1(t,phiS,t_ip);
    eta_ip  = interp1(t,eta ,t_ip);
    hf(i) = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('AHA ka=%.3g,M=%d,CPU=%.3g',ka,M,CPUTime));%[-1587 511 560 1000]
    if PLOT_WITH_SUBPLOT
        % subplot
        for iP = 1:nPannel
            subplot(nPannel,1,nPannel-iP+1), plot(x,eta_ip(iP,:),'k');hold on
            ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))))
            set(gca,'XTick',[],'XAxisLocation','origin');
            box off; grid on;
            % plot marker
            xDrift = mod(t_ip(iP)*c_p,L);
            yMarker =  interp1(x,eta_ip(iP,:),xDrift);
            hold on; plot(xDrift,yMarker,'ok','markersize',10)
            text(xDrift,yMarker,num2str(floor(t_ip(iP)*c_p/L)),'fontsize',10,'VerticalAlignment','middle','HorizontalAlignment','center')
        end
    else
        % single axes plot
        hold on;box off;
        set(gca,'XAxisLocation','origin');
        dz = 2.5*max(abs(eta_ip(1,:)));
        for iP = 1:nPannel
            z0 = (iP-1)*dz;
            plot(x,eta_ip(iP,:)+z0,'k','linewidth',1);
            text(L,z0,sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))));
            plot([0,L],[z0,z0],'k','linewidth',.5)
        end
        xDrift = mod(t*c_p,L);
        iBreak = find( diff(xDrift)<0 )+1; iBreak = [1;iBreak;length(t)];
        for iB=1:length(iBreak)-1, ii = iBreak([iB;iB+1])-[0;1];  plot(xDrift(ii),t(ii)/dt,'--r');end
        ylim([min(eta_ip(1,:)),2*max(eta_ip(1,:))+z0])
        xlabel('x [m]'); ylabel('\eta [m]');
    end
    drawnow
    if DO_EXPORT
       export_fig(hf(i),[filePath,cases{i},'_rePlot'],'-pdf'); 
    end
end