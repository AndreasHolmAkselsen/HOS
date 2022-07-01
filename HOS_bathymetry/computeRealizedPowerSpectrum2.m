% % script for making power spectra plots showing target spectrum and multiple simulations in overlay, all repeated in separate figures for
% varying sampling locations.

addpath c:/gits/timsas2/matlabLibs
clear

%% input
DO_EXPORT_S = 1; 

% matFiles={'./figures/mat/82100_closed_M5_H3p00_0p50_theta90_Nw1_nx2048_L50_Lb16p7_pad0_ikCut512_Md0_r0.mat'
%           './figures/mat/82100_closed_M5_H3p00_0p50_theta45_Nw1_nx2048_L50_Lb16p7_pad0_ikCut512_Md0_r0.mat'
%           './figures/mat/82100_closed_M5_H3p00_0p50_theta15_Nw1_nx2048_L50_Lb16p7_pad0_ikCut512_Md0_r0.mat'
% %           './figures/mat/82110_closed_M5_H0p50_theta_Nw1_nx2048_L50_Lb16p7_pad0_ikCut512_Md0_r0.mat'
% %             './figures/mat/82100_closed_M5_H3p00_0p50_theta90_Nw1_nx2048_L50_Lb16p7_pad0_ikCut1024_Md0_r0.mat'
% };
% legendStr = {'\theta = 90^o','\theta = 45^o','\theta = 15^o'};%,'Flat bottom'
% exportPath_S = './figures/powerSpec/82100/';


% matFiles={'./figures/mat/82100_closed_M5_H3p00_theta_Nw1_nx2048_L50_Lb16p7_pad0_ikCut1024_Md0_r0.mat'
%     './figures/mat/82110_closed_M5_H0p50_theta_Nw1_nx2048_L50_Lb16p7_pad0_ikCut512_Md0_r0.mat'};
% legendStr = {'3.0m water depth','0.5m water depth'};
% exportPath_S = './figures/powerSpec/82110_2xflatOnly/';

% matFiles={'./figures/mat/83000_closed_M5_H3p00_0p50_theta90_Nw1_nx2048_L50_Lb16p7_pad0_ikCut1024_Md0_r0.mat'
%           './figures/mat/83000_closed_M5_H3p00_0p50_theta45_Nw1_nx2048_L50_Lb16p7_pad0_ikCut1024_Md0_r0.mat'
%           './figures/mat/83000_closed_M5_H3p00_0p50_theta15_Nw1_nx2048_L50_Lb16p7_pad0_ikCut1024_Md0_r0.mat'
% };
% legendStr = {'\theta = 90^o','\theta = 45^o','\theta = 15^o'};
% exportPath_S = './figures/powerSpec/83000/';

matFiles={'./figures/mat/83000_closed_M5_H3p00_theta_Nw1_nx2048_L50_Lb16p7_pad0_ikCut1024_Md0_r0.mat'
          './figures/mat/83010_closed_M5_H0p50_theta_Nw1_nx2048_L50_Lb16p7_pad0_ikCut512_Md0_r0.mat'};
legendStr = {'3.0m water depth','0.5m water depth'};
exportPath_S = './figures/powerSpec/830x0_2xflatOnly/';

samplingLocations = 0:5:30;

%% code
nx_wp = length(samplingLocations);
% [hfS,haxS] = deal(zeros(nx_wp,1));
% plot first matFiles{1}
for iRN = 1:length(matFiles)
    load(matFiles{iRN})
    close(hf)
    % load(matFiles{1})
    for ix = 1:nx_wp, x_wp = samplingLocations(ix); % sampling location
        if iRN==1
            hfS(ix) = figure('color','w','name',['x = ',num2str(x_wp)]); 
            hfS(ix).Position(4) = 278;%hfS(ix).Position(3) = 1048;
            haxS(ix) = axes; grid on; hold on;
            waveDef = timsas_r1.data.WaveDefinition;
            waveDef.readFromSpecFile(waveMaker.signal{1}.specFile);
            STarg = waveDef.getSpectralDensity;
            plot(haxS(ix),waveDef.getUsedFrequencyVector,STarg,'k','linewidth',1.5,'DisplayName','Target');
            ylabel('S [m².s]');xlabel('frequency [Hz]');
%             title(rn);
            ylim([0,1.2*max(STarg)]);xlim([0,3.5/waveDef.component{1}.spectralDensity.Tp])
            legend('fontsize',10,'location','northeast'); %,'AutoUpdate','off'
        end
        
        % t_wp = t; % time is here uniform
        t_wp = waveMakerTime'; % sampling times
        dt = t_wp(2)-t_wp(1);
        % interpolate time signal h(x_wp,t)
        xx_wp = fzero(@(xx) real(fIP(xx*(xx>map.xxLR(1)&&xx<map.xxLR(2))))-x_wp,x_wp );
        if isnan(xx_wp), continue; end
        [~,i_wp] = min(abs(map.xi-xx_wp)); xx_wp = map.xi(i_wp); eta_wp = interp1(t,eta(:,i_wp),t_wp); % quick version
        %     [tt,xii] = meshgrid(t,map.xi);
        %     eta_wp = interp2(tt,xii,eta.',t_wp,xx_wp).';
        
        h_wp = imag( fIP(xx_wp+1i*eta_wp) ); % near enough in uniform regions.
        ts = timsas_r1.data.TimeSeries;
        ts.setFromVector(h_wp,'m',0,dt,'s');    
        ts.trim([50,inf]); % check!
        smoothingBandwidth = .01/waveDef.component{1}.spectralDensity.Tp;
        spec = ts.getSpectrum([0,99.9],1,smoothingBandwidth);
        plot(haxS(ix),spec.f,spec.S,'DisplayName',legendStr{iRN});%,'LineWidth',1.0);
    end
end


if DO_EXPORT_S
    if ~isfolder(exportPath_S),mkdir(exportPath_S);end
    for ix = 1:nx_wp, x_wp = samplingLocations(ix);
        fullPath = [exportPath_S,'/x_wp',num2str(x_wp)];
        savefig(hfS(ix),fullPath);
        export_fig(hfS(ix),fullPath,'-png','-pdf','-m2');
    end
end
