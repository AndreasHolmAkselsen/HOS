% name = 'waveField_Stokes_ka025_dt5T';
name = 'waveField_Stokes_ka028_dt2T';

hfAHA = open(['./doc/figures/AHA_',name,'.fig']);
hfSFo = open(['./doc/figures/SFo_',name,'.fig']);

hfSFo.Position=[527  0  1056  1000];

haAHA = hfAHA.Children;
haSFo = hfSFo.Children;

for i = 1:length(haSFo)
   hl = copyobj(haAHA(i).Children,haSFo(i));
   set(hl,'LineStyle','--','Color','r','LineWidth',1);
   haAHA(i).Children.LineWidth = 1;
end


savefig(hfSFo,['./doc/figures/AHAvsSFo_',name])
export_fig(hfSFo,['./doc/figures/AHAvsSFo_',name],'-pdf')