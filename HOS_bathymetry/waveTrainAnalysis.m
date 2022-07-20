% % script for making power spectra plots showing target spectrum and multiple simulations in overlay, all repeated in separate figures for
% varying sampling locations.

% addpath c:/gits/timsas2/matlabLibs
% clear
% close all
i=1;
%% input
DO_EXPORT_Aly = 1; 
exportPath_Aly = './figures/trainAly/';
fontSize = 11;

clear tStart nWaves nPeriods exPrefix refl matFiles

% 0.1
% tStart(i)=30;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_closed_M5_H1p00_0p10_theta90_NTdt5_nx4096_L152_Lb0_pad0_ikCut1024_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 0.2
% tStart(i)=30;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p20_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p1';i=i+1;
% % tStart(i)=35;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M10_H1p00_0p20_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p1';i=i+1;
% .3-.3 referenec
% tStart(i)=15;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H0p30_0p30_theta90_Nw3_nx2048_L76_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p32_ka0p2';i=i+1;
 %  0.3
% tStart(i)=30;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p2';i=i+1;
%       0.4, ka=.25 ( stable because of longer ramp and lower mode cut-off)
% tStart(i)=18;nWaves(i)=3;nPeriods(i)=3;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p40_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p25';i=i+1;
%       0.4, ka=.225 
% tStart(i)=18;nWaves(i)=3;nPeriods(i)=3;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p40_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p225';i=i+1;
% tStart(i)=18;nWaves(i)=3;nPeriods(i)=3;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p40_theta90_NTdt3_nx2048_L76_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p225';i=i+1;
%         0.5, ka=.25 
% tStart(i)=18;nWaves(i)=3;nPeriods(i)=3;exPrefix{i}='';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p50_theta90_NTdt3_nx2048_L76_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p25';i=i+1;



% 45 deg bathymetry, h2=0.1
% tStart(i)=30;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='45degDeep_';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p10_theta45_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 9 deg bathymetry, h2=0.1
% tStart(i)=35;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='9degDeep_';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p10_theta9_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 90deg to  h2=0.3 and then 9deg bathymetry, h3=0.1
% tStart(i)=35;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='9deg_';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_9_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 90deg to  h2=0.3 and then 4.5deg bathymetry, h3=0.1
% tStart(i)=35;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='4p5deg_';refl(i)=0;matFiles{i}='./figures/mat/train2_open_M5_H1p00_0p30_0p10_theta90_5_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 90deg to  h2=0.3 and then 2deg bathymetry, h3=0.1
% tStart(i)=40;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='2deg_';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_2_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 90deg to  h2=0.3 and then 1deg bathymetry, h3=0.1
% tStart(i)=40;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='1deg_';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_1_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;

% ka=0.075; 90deg to  h2=0.3 and then 1deg bathymetry, h3=0.1
% tStart(i)=40;nWaves(i)=10;nPeriods(i)=10;exPrefix{i}='1degka075_';refl(i)=0;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_1_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p075';i=i+1;



%%% REFLECTION WINDOWS
% 0.1
% tStart(i)=25;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl_';refl(i)=1;matFiles{i}='./figures/mat/train_closed_M5_H1p00_0p10_theta90_NTdt5_nx4096_L152_Lb0_pad0_ikCut1024_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% % tStart(i)=18;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl_';refl(i)=1;matFiles{i}='./figures/mat/train_closed_M5_H1p00_0p10_theta90_NTdt5_nx4096_L152_Lb0_pad0_ikCut1024_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 0.2
% tStart(i)=18;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p20_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p1';i=i+1;
% tStart(i)=18;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M10_H1p00_0p20_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p1';i=i+1;
 %  0.3
% tStart(i)=20;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p2';i=i+1;


% 45 deg bathymetry, h2=0.1
% tStart(i)=25;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl45degDeep_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p10_theta45_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% % tStart(i)=20;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl45degDeep_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p10_theta45_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 9 deg bathymetry, h2=0.1
% tStart(i)=25;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl9degDeep_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p10_theta9_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% % tStart(i)=20;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl9degDeep_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p10_theta9_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;

% 90deg to  h2=0.3 and then 9deg bathymetry, h3=0.1
% tStart(i)=40;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl9deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_9_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% % tStart(i)=35;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl9deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_9_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 90deg to  h2=0.3 and then 4.5deg bathymetry, h3=0.1
% tStart(i)=25;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl4p5deg_';refl(i)=1;matFiles{i}='./figures/mat/train2_open_M5_H1p00_0p30_0p10_theta90_5_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% % tStart(i)=20;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl4p5deg_';refl(i)=1;matFiles{i}='./figures/mat/train2_open_M5_H1p00_0p30_0p10_theta90_5_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 90deg to  h2=0.3 and then 1deg bathymetry, h3=0.1
% tStart(i)=25;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl2deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_2_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% % tStart(i)=20;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl2deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_2_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% 90deg to  h2=0.3 and then 1deg bathymetry, h3=0.1
% tStart(i)=25;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='refl1deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_1_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;
% % tStart(i)=20;nWaves(i)=8;nPeriods(i)=8;exPrefix{i}='refl1deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H1p00_0p30_0p10_theta90_1_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p05';i=i+1;

%%% inverse step aly, reflection
% 0.1
% % tStart(i)=17;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='downStrepRefl_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H0p10_1p00_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;
% tStart(i)=35;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='downStrepRefl_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H0p10_1p00_theta90_NTdt3p75_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;

% 0.2
% tStart(i)=32;nWaves(i)=4;nPeriods(i)=4;exPrefix{i}='downStrepRefl_';refl(i)=1;matFiles{i}='./figures/mat/train_closed_M5_H0p20_1p00_theta90_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T1p52_ka0p1';i=i+1;
% 0.3
% tStart(i)=20;nWaves(i)=4;nPeriods(i)=4;exPrefix{i}='downStrepRefl_';refl(i)=1;matFiles{i}='./figures/mat/train_closed_M5_H0p30_1p00_theta90_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T1p32_ka0p2';i=i+1;

% 0.1, 45 deg bathymetry
% tStart(i)=35;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='downStrepRefl45degDeep_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H0p10_1p00_theta45_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;
% 0.1, 9 deg bathymetry
tStart(i)=35;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='downStrepRefl9degDeep_';refl(i)=1;matFiles{i}='./figures/mat/train_closed_M5_H0p10_1p00_theta9_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;
% 0.1, 4.5 deg bathymetry
% tStart(i)=35;nWaves(i)=4;nPeriods(i)=4;exPrefix{i}='downStrepRefl4p5deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H0p10_0p30_1p00_theta90_5_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;
% 0.1, 2 deg bathymetry
% % tStart(i)=35;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='downStrepRefl2deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H0p10_0p30_1p00_theta90_2_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;
% tStart(i)=42;nWaves(i)=3;nPeriods(i)=3;exPrefix{i}='downStrepRefl2deg_';refl(i)=1;matFiles{i}='./figures/mat/train_open_M5_H0p10_0p30_1p00_theta90_2_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;
% 0.1, 1.0 deg bathymetry
% % tStart(i)=35;nWaves(i)=5;nPeriods(i)=5;exPrefix{i}='downStrepRefl1deg_';refl(i)=1;matFiles{i}='./figures/mat/train_closed_M5_H0p10_0p30_1p00_theta90_1_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;
% tStart(i)=42;nWaves(i)=3;nPeriods(i)=3;exPrefix{i}='downStrepRefl1deg_';refl(i)=1;matFiles{i}='./figures/mat/train_closed_M5_H0p10_0p30_1p00_theta90_1_NTdt3p63_nx4608_L150_Lb0_pad0_ikCut1152_Md0_r0_SSGW_T2p05_ka0p05';i=i+1;


%% code
for iRN = 1:length(matFiles)
    load(matFiles{iRN})
%     close(hf)
%     linkaxes(hf.Children)

    % use up-crossing on the first time point to determine wavelength/spatial grid
    
    iTStart = find(t>=tStart(iRN),1,'first');
    eta0 = eta(iTStart,:).';
    zS0 = fIP(map.xi+1i*eta0);
    xS0 = real(zS0); h = imag(zS0);
    dxStep = 0.0;
    crossShift = .5*max(h);
    if refl(iRN)
        xStep = min(real(zSingular(:)));
        iUC = find( h-crossShift<0&[h(2:end);0]-crossShift>0 & xS0<xStep-dxStep,nWaves(iRN)+1,'last' );
    else
        xStep = max(real(zSingular(:)));
        iUC = find( h-crossShift<0&[h(2:end);0]-crossShift>0 & xS0>xStep+dxStep,nWaves(iRN)+1,'first' );
    end
    iUC(1) = iUC(1)+1;
    x0 = xS0(iUC(1)-1) + (xS0(iUC(1))-xS0(iUC(1)-1))/(h(iUC(1))-h(iUC(1)-1))*(0-h(iUC(1)-1));
    xEnd = xS0(iUC(end)) + (xS0(iUC(end)+1)-xS0(iUC(end)))/(h(iUC(end)+1)-h(iUC(end)))*(0-h(iUC(end)));
    nxIP = 2*(iUC(end)-iUC(1)+1);
    dxIP = (xEnd-x0)/nxIP;
    xIP = x0+(0:nxIP-1)*dxIP;
    if refl(iRN)
        set(hf.Children,'XLim',[x0,xEnd+.75*(xEnd-x0)])
    else
        set(hf.Children,'XLim',[x0-.75*(xEnd-x0),xEnd])
    end
    
%     figure, plot(zS0(iUC(1)-1:iUC(end)+1));grid on;hold on;
%     plot([x0,xEnd],[0,0],'o')
    
    % interpolate data to uniform time based on up-crossing
    % (I have not bothered interpolating the time up-crossing point)
    ht = imag(fIP(map.xi(iUC(1)-1)+1i*eta(:,iUC(1)-1)));
%     ht = eta(:,iUC(1)-1);
    iUCt = find( ht-crossShift>0&[ht(2:end);0]-crossShift<0 & t>tStart(iRN),nPeriods(iRN)+1,'first' );
%     iUCt(1) = iUCt(1)+1;
    tEnd = t(iUCt(end));
    tStart_ = t(iUCt(1));
    ntIP = 1*(iUCt(end)-iUCt(1)+1);
    dtIP = (tEnd-tStart_)/ntIP;
    tIP = tStart_+(0:ntIP-1)'*dtIP;
    etaIP = interp1(t,eta,tIP);
    zSIP = fIP(map.xi.'+1i*etaIP);

    
    % interpolate sto uniform spatial grid on all times
    hIP = 0*(tIP+xIP);
    for i = 1:ntIP
        hIP(i,:) = interp1(real(zSIP(i,:)),imag(zSIP(i,:)),xIP);
%         hIP(i,:) = hIP(i,:)-mean(hIP(i,:));
    end    
    
    
    % fft vectors
    dw = 2*pi/(dtIP*ntIP);
    dk = 2*pi/(dxIP*nxIP);
    k0 = dk*nWaves(iRN);
    w0 = dw*nPeriods(iRN);
    K2 = findWaveNumbers(2*w0,map.H(end),0,0);
%     w = (0:ntIP-1)'*dw;
%     k = (0:nxIP-1)*dk;
%     fftt_h = fft(hIP,[],1);
%     fftx_h = fft(hIP,[],2);
%     ffttx_h = fft2(hIP);
    if mod(ntIP,2)==0
        w = [0:ntIP/2-1,-ntIP/2:-1]'*dw;
    else
        w = [0:(ntIP-1)/2, -(ntIP-1)/2:-1]'*dw;   
    end
    if mod(nxIP,2)==0
        k = [0:nxIP/2-1,-nxIP/2:-1]'*dk;
    else
        k = [0:(nxIP-1)/2, -(nxIP-1)/2:-1]'*dk;   
    end
    w = fftshift(w); k = fftshift(k);
    fftt_h = fftshift(fft(hIP,[],1),1);
    fftx_h = fftshift(fft(hIP,[],2),2);
    ffttx_h = fftshift(fft2(hIP));
        
    
    hfAly.xt = figure('color','w'); 
%     [xx,tt]=meshgrid(xIP,tIP);surf(xx,tt,hIP,'LineStyle','none');xlabel('x');ylabel('t')
    imagesc(xIP,tIP,hIP);xlabel('Position','FontSize',fontSize);ylabel('Time','FontSize',fontSize), axis xy
    
    
    hfAly.kt = figure('color','w'); 
    imagesc(k/k0,tIP,abs(fftx_h));xlabel('k/k_0','FontSize',fontSize);ylabel('Time','FontSize',fontSize), axis xy

    hfAly.xw = figure('color','w'); 
    imagesc(xIP,w/w0,abs(fftt_h));xlabel('Position','FontSize',fontSize);ylabel('\omega/\omega_0','FontSize',fontSize), axis xy
    
    hfAly.kw = figure('color','w'); 
    ffttx_h(1)=nan;ii = 20;
%     imagesc(k(1:ii)/k0,w(1:ii)/w0,abs(ffttx_h(1:ii,1:ii)));xlabel('k/k_0');ylabel('\omega/\omega_0'), axis xy
    imagesc(k/k0,w/w0,abs(ffttx_h));xlabel('k/k_0','FontSize',fontSize);ylabel('\omega/\omega_0','FontSize',fontSize), axis xy
    
    if refl(iRN) 
%         xlim(hfAly.kt.Children,[0,3.5])
%         ylim(hfAly.xw.Children,[0,3.5])    
%         axis(hfAly.kw.Children,[-3.5,3.5,-2.5,2.5])
        xlim(hfAly.kt.Children,[0,5.5])
        ylim(hfAly.xw.Children,[0,5.5])    
        axis(hfAly.kw.Children,[-5.5,5.5,-5.5,5.5])
    else
        xlim(hfAly.kt.Children,[0,5])
        ylim(hfAly.xw.Children,[0,5])    
        axis(hfAly.kw.Children,[0,5,-4,0])
    end
    
    
%     figure; 
%     phaseAngle = unwrap(angle(fftx_h(:,nWaves+1)));
% %     plot(tIP,phaseAngle);
%     plot(.5*(tIP(1:end-1)+tIP(2:end)),-diff(phaseAngle)/dtIP);hold on; grid on
%     plot(tIP([1,end]),w0*[1,1],'k',tIP([1,end]),2*pi/IC.T*[1,1],'k')
  

    if DO_EXPORT_Aly
        if ~isfolder(exportPath_Aly),mkdir(exportPath_Aly);end
        figTypes = {'xt','kt','xw','kw'};
        Hstr = sprintf('%.1f_',map.H);
        fileName = sprintf('%sH%ska%.3g_nw%d_np%d',exPrefix{iRN},Hstr,IC.ka,nWaves(iRN),nPeriods(iRN));
        fileName(fileName=='.') = 'p';
        for i=1:length(figTypes)
            ft = figTypes{i};
            fullPath = [exportPath_Aly,'/',fileName,'_',ft];
            tic
            savefig(hfAly.(ft),fullPath);
            fprintf('Time %s, savefig: %.3gs\n',ft,toc)
            tic
            export_fig(hfAly.(ft),fullPath,'-png','-pdf','-m2');
            fprintf('Time %s, export_fig: %.3gs\n',ft,toc)
        end
        tic
        export_fig(hf,[exportPath_Aly,'/',fileName],'-png','-pdf','-m2');
        fprintf('Time resave: %.3gs\n',toc)
    end

end


