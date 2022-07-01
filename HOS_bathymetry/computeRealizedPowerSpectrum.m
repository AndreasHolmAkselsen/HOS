
% load('./figures/mat/82100_closed_M5_H3p00_0p50_theta90_Nw1_nx2048_L50_Lb16p7_pad0_ikCut512_Md0_r0.mat')


hfS = figure('color','w');axes;
ylabel(['S [',spec.spectrumUnit.getUnitName,']'])
xlabel('frequency [Hz]')
title([rn]);
ylim([0,8e-3]);xlim([0,1.8])
grid on; hold on;
ts = timsas_r1.data.TimeSeries;
ts.setFromVector(h_wp,'m',0,dt,'s');


[tt,xii] = meshgrid(t,map.xi);

% t_wp = t; % time is here uniform
t_wp = waveMakerTime'; % sampling times
dt = t_wp(2)-t_wp(1);
% smoothingBandwidth = 0;
% [f,S] = level1.spec.computePowerSpectralDensity(eta_wp,dt,smoothingBandwidth);
waveDef = timsas_r1.data.WaveDefinition;
waveDef.readFromSpecFile(waveMaker.signal{1}.specFile);
plot(waveDef.getUsedFrequencyVector,waveDef.getSpectralDensity,'r','linewidth',1.5,'DisplayName','Target');

for x_wp = 0:5:30 % sampling location
    
    % interpolate time signal h(x_wp,t)
    xx_wp = fzero(@(xx) real(fIP(xx*(xx>map.xxLR(1)&&xx<map.xxLR(2))))-x_wp,x_wp );
    if isnan(xx_wp), continue; end
%     [~,i_wp] = min(abs(map.xi-xx_wp)); xx_wp = map.xi(i_wp); eta_wp = interp1(t,eta(:,i_wp),t_wp); % quick version
    eta_wp = interp2(tt,xii,eta.',t_wp,xx_wp).';

    h_wp = imag( fIP(xx_wp+1i*eta_wp) ); % near enough in uniform regions.    


    ts.value = h_wp;
    ts.trim([50,inf]); % check!
    
    smoothingBandwidth = .01/waveDef.component{1}.spectralDensity.Tp;
    spec = ts.getSpectrum([0,99.9],1,smoothingBandwidth);
    plot(spec.f,spec.S,'DisplayName',['x = ',num2str(x_wp)]);
end
legend('fontsize',10,'location','northeast');

% if DO_EXPORT
%     savefig(hfS,[exportPath,'/fig/S_',fileName]);
%     export_fig(hfS,[exportPath,'/S_',num2str(x_wp,'%.3g'),'_',fileName],exportFormats{:});
% end


% x_wp = 0; % sampling location
% t_wp = param.waveMaker.time'; % sampling times
% 
% % interpolate time signal h(x_wp,t) 
% xx_wp_temp = fzero(@(xx) real(fIP(xx))-x_wp,x_wp );
% [~,i_wp] = min(abs(map.xi-xx_wp_temp)); 
% xx_wp = map.xi(i_wp);
% eta_wp = interp1(t,eta(:,i_wp),t_wp); 
% h_wp = imag( fIP(xx_wp+1i*eta_wp) ); % near enough in uniform regions.
% 
% % [~,i_wp] = min(abs(x-x_wp)); 
% % h_wp = interp1(t,eta(:,i_wp),t_uniform); 
% 
% 
% dt = t_wp(2)-t_wp(1);
% % smoothingBandwidth = 0;
% % [f,S] = level1.spec.computePowerSpectralDensity(eta_wp,dt,smoothingBandwidth);
% waveDef = timsas_r1.data.WaveDefinition;
% waveDef.readFromSpecFile(waveMaker.signal{1}.specFile);
% ts = timsas_r1.data.TimeSeries;
% ts.setFromVector(h_wp,'m',0,dt,'s');
% ts.plot
% ts.trim([50,inf]);
% 
% smoothingBandwidth = .01/waveDef.component{1}.spectralDensity.Tp;
% spec = ts.getSpectrum([1,99.9],1,smoothingBandwidth);
% 
% figure, plot(spec.f,spec.S,waveDef.getUsedFrequencyVector,waveDef.getSpectralDensity);grid on;
% ylabel(['S [',spec.spectrumUnit.getUnitName,']'])
% xlabel('frequency [Hz]')
