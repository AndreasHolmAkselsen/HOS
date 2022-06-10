x_wp = 0; % sampling location
t_wp = param.waveMaker.time'; % sampling times

% interpolate time signal h(x_wp,t) 
xx_wp_temp = fzero(@(xx) real(fIP(xx))-x_wp,x_wp );
[~,i_wp] = min(abs(map.xi-xx_wp_temp)); 
xx_wp = map.xi(i_wp);
eta_wp = interp1(t,eta(:,i_wp),t_wp); 
h_wp = imag( fIP(xx_wp+1i*eta_wp) ); % near enough in uniform regions.

% [~,i_wp] = min(abs(x-x_wp)); 
% h_wp = interp1(t,eta(:,i_wp),t_uniform); 


dt = t_wp(2)-t_wp(1);
% smoothingBandwidth = 0;
% [f,S] = level1.spec.computePowerSpectralDensity(eta_wp,dt,smoothingBandwidth);
waveDef = timsas_r1.data.WaveDefinition;
waveDef.readFromSpecFile(waveMaker.signal{1}.specFile);
ts = timsas_r1.data.TimeSeries;
ts.setFromVector(h_wp,'m',0,dt,'s');
ts.plot
ts.trim([50,inf]);

smoothingBandwidth = .01/waveDef.component{1}.spectralDensity.Tp;
spec = ts.getSpectrum([1,99.9],1,smoothingBandwidth);

figure, plot(spec.f,spec.S,waveDef.getUsedFrequencyVector,waveDef.getSpectralDensity);grid on;
ylabel(['S [',spec.spectrumUnit.getUnitName,']'])
xlabel('frequency [Hz]')
