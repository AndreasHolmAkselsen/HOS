% % script for making power spectra plots showing target spectrum and multiple simulations in overlay, all repeated in separate figures for
% varying sampling locations.

addpath c:/gits/timsas2/matlabLibs
% clear
close all

%% input
DO_EXPORT_Aly = 0; 
exportPath_Aly = './figures/trainAly/';
matFiles={
    % 0.1
%     './figures/mat/train_closed_M5_H1p00_0p10_theta90_NTdt5_nx4096_L152_Lb0_pad0_ikCut1024_Md0_r0_SSGW_T1p13_ka0p05'
    % 0.2
%     './figures/mat/train_open_M5_H1p00_0p20_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut768_Md0_r0_SSGW_T1p13_ka0p1'
        % .3-.3 referenec
%           './figures/mat/train_open_M5_H0p30_0p30_theta90_Nw3_nx2048_L76_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p32_ka0p2'
        %  0.3
%           './figures/mat/train_open_M5_H1p00_0p30_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p2'
%       0.4, ka=.25 ( stable because of longer ramp and lower mode cut-off)
%             './figures/mat/train_open_M5_H1p00_0p40_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p25'
          %       0.4, ka=.225 
%           './figures/mat/train_open_M5_H1p00_0p40_theta90_NTdt5_nx3072_L100_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p225'
%           './figures/mat/train_open_M5_H1p00_0p40_theta90_NTdt3_nx2048_L76_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p225'
%         0.5, ka=.25 
          './figures/mat/train_open_M5_H1p00_0p50_theta90_NTdt3_nx2048_L76_Lb0_pad0_ikCut512_Md0_r0_SSGW_T1p13_ka0p25'

        };



tStart = 18; nWaves = 3; nPeriods = 3;
% tStart = 15; nWaves = 5; nPeriods = 5;
% tStart = 30; nWaves = 10; nPeriods = 10;
% tStart = 30; nWaves = 8; nPeriods = 8;
% tStart = 25; nWaves = 5; nPeriods = 5;
%% code
for iRN = 1:length(matFiles)
    load(matFiles{iRN})
%     close(hf)
%     linkaxes(hf.Children)

    % use up-crossing on the first time point to determine wavelength/spatial grid
    
    iTStart = find(t>=tStart,1,'first');
    tStart = t(iTStart);

    eta0 = eta(iTStart,:).';
    zS0 = fIP(map.xi+1i*eta0);
    h = imag(zS0);
    xStep = max(real(zSingular));
    dxStep = 0.0;
    iUC = find( h<0&[h(2:end);0]>0 & real(zS0)>xStep+dxStep,nWaves+1,'first' );
    iUC(1) = iUC(1)+1;
    x0 = real(zS0(iUC(1)-1)) + real(zS0(iUC(1))-zS0(iUC(1)-1))/imag(zS0(iUC(1))-zS0(iUC(1)-1))*(0-imag(zS0(iUC(1)-1)));
    xEnd = real(zS0(iUC(end))) + real(zS0(iUC(end)+1)-zS0(iUC(end)))/imag(zS0(iUC(end)+1)-zS0(iUC(end)))*(0-imag(zS0(iUC(end))));
    nxIP = 2*(iUC(end)-iUC(1)+1);
    dxIP = (xEnd-x0)/nxIP;
    xIP = x0+(0:nxIP-1)*dxIP;
    set(hf.Children,'XLim',[x0-.5*(xEnd-x0),xEnd])
%     figure, plot(zS0(iUC(1)-1:iUC(end)+1));grid on;hold on;
%     plot([x0,xEnd],[0,0],'o')
    
    % interpolate data to uniform time based on up-crossing
    % (I have not bothered interpolating the time up-crossing point)
    ht = imag(fIP(map.xi(iUC(1)-1)+1i*eta(:,iUC(1)-1)));
%     ht = eta(:,iUC(1)-1);
    iUCt = find( ht>0&[ht(2:end);0]<0 & t>tStart,nPeriods+1,'first' );
%     iUCt(1) = iUCt(1)+1;
    tEnd = t(iUCt(end));
    tStart = t(iUCt(1));
    ntIP = 1*(iUCt(end)-iUCt(1)+1);
    dtIP = (tEnd-tStart)/ntIP;
    tIP = tStart+(0:ntIP-1)'*dtIP;
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
    k0 = dk*nWaves;
    w0 = dw*nPeriods;
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
    imagesc(xIP,tIP,hIP);xlabel('Position');ylabel('Time'), axis xy
    
    
    hfAly.kt = figure('color','w'); 
    imagesc(k/k0,tIP,abs(fftx_h));xlabel('k/k_0');ylabel('Time'), axis xy
    xlim([0,5])

    hfAly.xw = figure('color','w'); 
    imagesc(xIP,w/w0,abs(fftt_h));xlabel('Position');ylabel('\omega/\omega_0'), axis xy
    ylim([0,5])
    
    hfAly.kw = figure('color','w'); 
    ffttx_h(1)=nan;ii = 20;
%     imagesc(k(1:ii)/k0,w(1:ii)/w0,abs(ffttx_h(1:ii,1:ii)));xlabel('k/k_0');ylabel('\omega/\omega_0'), axis xy
    imagesc(k/k0,w/w0,abs(ffttx_h));xlabel('k/k_0');ylabel('\omega/\omega_0'), axis xy
    axis([0,5,-4,0])
    
%     figure; 
%     phaseAngle = unwrap(angle(fftx_h(:,nWaves+1)));
% %     plot(tIP,phaseAngle);
%     plot(.5*(tIP(1:end-1)+tIP(2:end)),-diff(phaseAngle)/dtIP);hold on; grid on
%     plot(tIP([1,end]),w0*[1,1],'k',tIP([1,end]),2*pi/IC.T*[1,1],'k')
  

    if DO_EXPORT_Aly
        if ~isfolder(exportPath_Aly),mkdir(exportPath_Aly);end
        figTypes = {'xt','kt','xw','kw'};
        fileName = sprintf('H%.1f_%.1f_ka%.3g_nw%d_np%d',map.H(1),map.H(2),IC.ka,nWaves,nPeriods);
        fileName(fileName=='.') = 'p';
        for i=1:length(figTypes)
            ft = figTypes{i};
            fullPath = [exportPath_Aly,'/',fileName,'_',ft];
            savefig(hfAly.(ft),fullPath);
            export_fig(hfAly.(ft),fullPath,'-png','-pdf','-m2');
        end
        export_fig(hf,[exportPath_Aly,'/',fileName],'-png','-pdf','-m2');
    end

end


