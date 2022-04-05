
hf(1) = openfig(['./figures/testConf__ka0p05_M3_H1p00_0p99_Nw30_dt2p35T_nx1920_pad0_kCut400.fig']);
hf(2) = openfig(['./figures/testFlat_Taylor_ka0p05_M3_h1p00_Nw30_dt2p35T_nx1920_pad0_kCut400.fig'],'invisible');
legendNames = {'Conformal','flat'};


lineStyles = {'-','--',':'};
colors = {'k','r','b'};
hf(1).Position=[527  0  1056  1000];

ha0 = hf(1).Children;
for iF = 2:numel(hf)
    ha = hf(iF).Children;
    for i = 1:length(ha0)
        hl = copyobj(ha(i).Children,ha0(i));
        set(hl,'LineStyle',lineStyles{iF},'Color',colors{iF},'LineWidth',1);
    end
    close(hf(iF));
end
legend(ha0(end),legendNames)


return

exportPath = ['./figures/vortex'];
savefig(hf(1),exportPath)
export_fig(hf(1),exportPath,'-pdf','-png')