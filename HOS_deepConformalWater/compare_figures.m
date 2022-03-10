% name = 'waveField_Stokes_ka025_dt5T';
name = 'ka0p25_M5_Nw10_dt1T_nx512';

% prefix = {'Taylor','phiS_x','phi_x'};
% for i = 1:length(prefix)
%     hf(i) = openfig(['./figures/',prefix{i},'ka0p28_M1_Nw10_dt1T_nx512.fig'],'invisible');
% end
% hf(1).Visible = 'on';

hf(1) = openfig(['./figures/Taylorka0p28_M5_Nw10_dt1T_nx512.fig']);
hf(2) = openfig(['./figures/phiS_xka0p28_M1_Nw10_dt1T_nx512.fig'],'invisible');
hf(3) = openfig(['./figures/phi_xka0p28_M1_Nw10_dt1T_nx512.fig'],'invisible');
legendNames = {'Taylor','phiS_x','phi_x'};
lineStyles = {'-','--',':'};
colors = {'k','r','b'}
hf(1).Position=[527  0  1056  1000];

ha0 = hf(1).Children;
for iF = 2:numel(hf)
    ha = hf(iF).Children;
    for i = 1:length(ha0)
        hl = copyobj(ha(i).Children,ha0(i));
        set(hl,'LineStyle',lineStyles{iF},'Color',colors{iF},'LineWidth',1);
        ha(i).Children.LineWidth = 1;
    end
    close(hf(iF));
end
legend(ha0(1),legendNames)
% 
% haAHA = hfAHA.Children;
% haSFo = hfSFo.Children;
% 
% for i = 1:length(haSFo)
%    hl = copyobj(haAHA(i).Children,haSFo(i));
%    set(hl,'LineStyle','--','Color','r','LineWidth',1);
%    haAHA(i).Children.LineWidth = 1;
% end
% close(hfAHA)
% legend(haSFo(1),{'SFo','AHA'})

% savefig(hfSFo,['./HOS_SFo_curr/figures/AHAvsSFo_',name])
% export_fig(hfSFo,['./HOS_SFo_curr/figures/AHAvsSFo_',name],'-pdf')