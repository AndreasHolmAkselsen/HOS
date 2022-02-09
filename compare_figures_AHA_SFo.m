% name = 'waveField_Stokes_ka025_dt5T';
name = 'ka0p25_M5_Nw10_dt1T_nx512';

hfSFo = openfig(['./HOS_SFo_curr\figures/SFo_',name,'.fig']);
hfAHA = openfig(['./HOS_SFo_curr/figures/AHA_',name,'.fig'],'invisible');

hfSFo.Position=[527  0  1056  1000];

haAHA = hfAHA.Children;
haSFo = hfSFo.Children;

for i = 1:length(haSFo)
   hl = copyobj(haAHA(i).Children,haSFo(i));
   set(hl,'LineStyle','--','Color','r','LineWidth',1);
   haAHA(i).Children.LineWidth = 1;
end
close(hfAHA)
legend(haSFo(1),{'SFo','AHA'})

savefig(hfSFo,['./HOS_SFo_curr/figures/AHAvsSFo_',name])
export_fig(hfSFo,['./HOS_SFo_curr/figures/AHAvsSFo_',name],'-pdf')