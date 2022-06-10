

x_wp = 0;
[~,i_wp] = min(abs(x-x_wp)); 

t_uniform = param.waveMaker.time';
eta_wp = interp1(t,eta(:,i_wp),t_uniform); 


dt = t_uniform(2)-t_uniform(1);
% smoothingBandwidth = 0;
% [f,S] = level1.spec.computePowerSpectralDensity(eta_wp,dt,smoothingBandwidth);

waveDef = timsas_r1.data.WaveDefinition;
waveDef.readFromSpecFile(waveMaker.signal{1}.specFile);

ts = timsas_r1.data.TimeSeries;
ts.setFromVector(eta_wp,'m',0,dt,'s');

ts.plot
ts.trim([50,inf]);

smoothingBandwidth = .01/waveDef.component{1}.spectralDensity.Tp;
spec = ts.getSpectrum([1,99.9],1,smoothingBandwidth);

figure, plot(spec.f,spec.S,waveDef.getUsedFrequencyVector,waveDef.getSpectralDensity);grid on;
ylabel(['S [',spec.spectrumUnit.getUnitName,']'])
xlabel('frequency [Hz]')