clear variables
close all
cCol={'b','r','g','k','m'};
xWind=[-6 6];
% tWind=[-1 -0.5 0 0.5 1];
tWind=0.4;

%Path to wave data
dataPath='.\..\simResults';

%Simulations to compare
simList='testHOS_BichromWave_O4';

%Main loop
f=figure('Position',[1500,50,700,200]);
sp=cell(1,length(tWind));

%Load data
simResFileName=[dataPath filesep simList '_simNWT.mat'];
temp=load(simResFileName);simRes=temp.simRes;
velResFileName=[dataPath filesep simList '_velNWT.mat'];
velRes=load(velResFileName);

%Find time and location of steepest wave
xslMax=77.2283;
tslMax=69.6600;
xmin=xslMax+xWind(1);
xmax=xslMax+xWind(2);

k1=simRes.nwtSpec.init.waveComp{1}.nk0*2*pi/simRes.nwtSpec.sim.Lx;
k2=simRes.nwtSpec.init.waveComp{2}.nk0*2*pi/simRes.nwtSpec.sim.Lx;
km=(k1+k2)/2;
sig1=sqrt(9.81*k1);
sig2=sqrt(9.81*k2);
sigm=(sig1+sig2)/2;
Tm=2*pi/sigm;
for iSP=1:length(tWind)
    [~,indtVel]=min(abs(velRes.tUsed-(tslMax+tWind(iSP))));
    [~,indt]=min(abs(simRes.t-(velRes.tUsed(indtVel))));
    sp{iSP}=subplot(length(tWind),1,iSP);
    set(sp{iSP},'nextplot','replacechildren');
    imagesc(km*velRes.xOut,km*velRes.zOut,km/sigm*velRes.u(:,:,indtVel).','AlphaData',~isnan(velRes.u(:,:,indtVel).'))
    hold on
    plot(km*simRes.x,km*simRes.eta(indt,:),'k')

%     imagesc(km*velRes.xOut,km*velRes.zOut,km/sigm*velRes.u(:,:,indtVel).','AlphaData',~isnan(velRes.u(:,:,indtVel).'))
%     hold on
%     plot(km*(simRes.alpha+simRes.x(indt,:)),km*simRes.z(indt,:),'k','linewidth',1)
    %         set(gca,'XLim',[xmin xmax])
end


%Finalize plot
for iSP=1:length(tWind)
    subplot(length(tWind),1,iSP)
    grid on
    axis equal
    title(['$t/T_m=' num2str((tslMax+tWind(iSP))/Tm,'%6.2f') '$'],'interpreter','latex')
    ylabel('$k_m \hspace{0.05cm} z$','interpreter','latex')
    if iSP==length(tWind)
        xlabel('$k_m \hspace{0.05cm} x$','interpreter','latex');
    end
    set(gca,'FontSize',9)
    %     caxis(km/sigm*[-1 1.6])
    caxis([-0.3 0.6])
    cb=colorbar;
    if iSP==1
        set(get(cb,'title'),'string','$\frac{k_m}{\sigma_m} u$','interpreter','latex');
    end
    set(get(cb,'title'),'FontSize',12);    colormap jet
    set(sp{iSP},'XLim',[km*xmin km*xmax])
    ylim([-1*km 0.5])
    if iSP~=length(tWind)
        set(sp{iSP},'XTickLabel','')
    end
end
saveas(f,'.\..\rapPlots\VelHOS_O4','png')
% set(gcf,'renderer','painters')
set(gcf,'renderer','opengl')
print(f,'.\..\..\..\Paper\ComputationalPhysics\melhos\VelHOS_O4','-depsc')

