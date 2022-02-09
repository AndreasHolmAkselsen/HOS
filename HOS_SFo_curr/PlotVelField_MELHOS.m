clear variables
close all
cCol={'b','r','g','k','m'};
xWind=[-6 6];
% tWind=[-1 -0.5 0 0.5 1];
tWind=-1.2:0.4:1.2;

%Path to wave data
dataPath='.\..\simResults';

%Simulations to compare
simList={'testMELHOS_BichromWave_O4'};

%Main loop
f=figure('Position',[1500,50,700,800]);
% for iSP=1:length(tWind)
%     subplot(length(tWind),1,iSP)
%     hold on
% end
sp=cell(1,length(tWind));
for iSim=1:length(simList)
    
    %Load data
    simResFileName=[dataPath filesep simList{iSim} '_simNWT.mat'];
    temp=load(simResFileName);simRes=temp.simRes;
    velResFileName=[dataPath filesep simList{iSim} '_velNWT.mat'];
    velRes=load(velResFileName);
    
    %Find time and location of steepest wave
    [alphaMat,tMat]=meshgrid(simRes.alpha,simRes.t);
    dalpha=simRes.alpha(2)-simRes.alpha(1);
    slope=([simRes.z(:,2:end) simRes.z(:,1)]-[simRes.z(:,end) simRes.z(:,1:end-1)])./...
        (2*dalpha+[simRes.x(:,2:end) simRes.x(:,1)]-[simRes.x(:,end) simRes.x(:,1:end-1)]);
    [slmax,indslMax]=max(-slope(:));
    if iSim==1
        indtslMax=find(simRes.t==tMat(indslMax));
        indaslMax=find(simRes.alpha==alphaMat(indslMax));
        xslMax=simRes.alpha(indaslMax)+simRes.x(indtslMax,indaslMax);
        tslMax=simRes.t(indtslMax);
        xmin=xslMax+xWind(1);
        xmax=xslMax+xWind(2);
        k1=simRes.nwtSpec.init.waveComp{1}.nk0*2*pi/simRes.nwtSpec.sim.Lx;
        k2=simRes.nwtSpec.init.waveComp{2}.nk0*2*pi/simRes.nwtSpec.sim.Lx;
        km=(k1+k2)/2;
        sig1=sqrt(9.81*k1);
        sig2=sqrt(9.81*k2);
        sigm=(sig1+sig2)/2;
        Tm=2*pi/sigm;
        Lx=simRes.nwtSpec.sim.Lx;
    end
    for iSP=1:length(tWind)
        [~,indtVel]=min(abs(velRes.tUsed-(tslMax+tWind(iSP))));
        [~,indt]=min(abs(simRes.t-(velRes.tUsed(indtVel))));
        sp{iSP}=subplot(length(tWind),1,iSP);
        set(sp{iSP},'nextplot','replacechildren');
        imagesc(km*velRes.xOut,km*velRes.zOut,km/sigm*velRes.u(:,:,indtVel).','AlphaData',~isnan(velRes.u(:,:,indtVel).'))
        hold on
        plot(km*(simRes.alpha+simRes.x(indt,:)),km*simRes.z(indt,:),'k','linewidth',1)
%         set(gca,'XLim',[xmin xmax])
    end
    
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
saveas(f,'.\..\rapPlots\VelMELHOS_O3','png')
% set(gcf,'renderer','painters')
set(gcf,'renderer','opengl')
print(f,'.\..\..\..\Paper\ComputationalPhysics\melhos\VelHOSMELHOS_O3','-depsc')

%Plot time series of surface elevation at the center of the domain
f2=figure('Position',[1600 700 900 150]);
plot(simRes.t/Tm,simRes.eta(:,798)*km,'b')
grid on
xlim([0 80])
xlabel('$t/T_m$','interpreter','latex')
ylabel('$k_m \hspace{0.05cm} \eta$','interpreter','latex')
print(f2,'.\..\..\..\Paper\ComputationalPhysics\melhos\SurfElev_O4','-depsc')

