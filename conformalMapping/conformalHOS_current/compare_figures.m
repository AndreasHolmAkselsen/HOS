
% prefix = {'Taylor','phiS_x','phi_x'};
% for i = 1:length(prefix)
%     hf(i) = openfig(['./figures/',prefix{i},'ka0p28_M1_Nw10_dt1T_nx512.fig'],'invisible');
% end
% hf(1).Visible = 'on';
% 

hf(1) = openfig(['./figures/Taylorka0p4_M5_h100p00_Nw10_dt2p5T_nx1024_pad0_ikCut5.fig']);
hf(2) = openfig(['./figures/Taylorka0p4_M5_h100p00_Nw10_dt2p5T_nx1024_pad1_ikCut5.fig'],'invisible');
hf(3) = openfig(['./figures/SFo_ka0p4_M5_h100p00_nx1024_ikCut50.fig'],'invisible');
legendNames = {'unpadded','padded','SFo'};


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