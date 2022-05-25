
hf(1) = openfig(['./figures/fig/closed2_zeroMap_SSGW_T2p50_ka0p1_M5_H1p00_theta_Nw30_dt2T_nx1920_pad0_ikCut300_Md0_r0.fig']);
hf(2) = openfig(['./figures/fig/closed_zeroMap_SSGW_T2p50_ka0p1_M5_H1p00_theta_Nw30_dt2T_nx1920_pad0_ikCut300_Md0_r0.fig'],'invisible');
% hf(2) = openfig(['../../HOS_ChalikovTaylor_current/figures/benchmarkBathimetry_Taylor_ka0p1_M5_h1p00_Nw30_dt2T_nx1920_pad0_kCut8p993.fig'],'invisible');
hf(3) = openfig(['./figures/fig/closed3_zeroMap_SSGW_T2p50_ka0p1_M5_H1p00_theta_Nw30_dt2T_nx1920_pad0_ikCut300_Md0_r0.fig'],'invisible');

legendNames = {'closed','new L','new L, alternative fft'};


lineStyles = {'-','--',':'};
colors = {'k','r','b'};
hf(1).Position=[527  0  1056  1000];

ha0 = hf(1).Children;
for iF = 2:numel(hf)
    ha = hf(iF).Children;
    for i = 1:length(ha0)
        hl = copyobj(ha(i).Children,ha0(i));
%         hl.XData = hl.XData-mean(hl.XData);
        set(hl,'LineStyle',lineStyles{iF},'Color',colors{iF},'LineWidth',1);
    end
    close(hf(iF));
end
legend(ha0(end),legendNames)

hf2 = figure('color','w');
ha2=copyobj(ha0(end),hf2); ha2.Position([2,4])=[.15,.8];

return

exportPath = ['./figures/vortex'];
savefig(hf(1),exportPath)
export_fig(hf(1),exportPath,'-pdf','-png')