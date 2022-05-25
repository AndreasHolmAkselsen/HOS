clear hf;
hf(1) = openfig(['./figures/fig/closed2_zeroMap_SSGW_T2p50_ka0p1_M5_H1p00_theta_Nw30_dt2T_nx1920_pad0_ikCut300_Md0_r0.fig']);
hf(2) = openfig(['./figures/fig/closed_zeroMap_closed_SSGW_T2p50_ka0p1_M5_H1p00_theta_Nw30_dt2T_nx1920_pad0_ikCut300_Md0_r0.fig'],'invisible');
hf(3) = openfig(['./figures/fig/closed_zeroMap_open_SSGW_T2p50_ka0p1_M5_H1p00_theta_Nw30_dt2T_nx1920_pad0_ikCut300_Md0_r0.fig'],'invisible');
hf(4) = openfig(['../../figures/fig/open_zeroMap_SSGW_T2p50_ka0p1_M5_H1p00_theta_Nw30_dt2T_nx1920_pad0_ikCut300_Md0_r0.fig'],'invisible');

legendNames = {'before merger','after merger','open domain after merger','old open domain'};

lineStyles = {'-','--',':','-.'};
colors = {'k','r','b','c'};
hf(1).Position=[527  0  1056  1000];


ha = hf(1).Children;
ha0 = hf(1).Children;
for i = 1:length(ha0)
    for iChild = 1:length(ha(i).Children)
        if length(ha(i).Children(iChild).XData) <= 2, continue; end
        hl = ha(i).Children(iChild);break
    end
    set(hl,'LineStyle',lineStyles{1},'Color',colors{1},'LineWidth',1,'displayName',legendNames{1});
end


for iF = 2:numel(hf)
    ha = hf(iF).Children;
    for i = 1:length(ha0)
        for iChild = 1:length(ha(i).Children)
            if length(ha(i).Children(iChild).XData) <= 2, continue; end
            hl = copyobj(ha(i).Children(iChild),ha0(i));break
        end
%         hl.XData = hl.XData-mean(hl.XData);
        set(hl,'LineStyle',lineStyles{iF},'Color',colors{iF},'LineWidth',1,'displayName',legendNames{iF});
    end
    close(hf(iF));
end
legend(ha0(end))

hf2 = figure('color','w');
ha2=copyobj(ha0(end),hf2);  ha2.Position([2,4])=[.15,.8];

return

exportPath = ['./figures/'];
savefig(hf(1),exportPath)
export_fig(hf(1),exportPath,'-pdf','-png')