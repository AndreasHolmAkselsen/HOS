clear variables 
close all
cCol={'b','r','g','k','m'};
xWind=[-6 6];
% tWind=[-1 -0.5 0 0.5 1];
tWind=-1.2:0.4:1.2;

%Path to wave data
dataPath='.\..\simResults';

%Simulations to compare
% simList={'testMELHOS_BichromWave','testHOS_BichromWave'};
simList={'testMELHOS_BichromWave_O4','testHOS_BichromWave_O4'};

%Main loop
f=figure('Position',[1500,50,1200,600]);
for iSP=1:length(tWind)
    subplot(length(tWind),1,iSP)
    hold on
end
for iSim=[1:length(simList),1]

    %Load data
    simResFileName=[dataPath filesep simList{iSim} '_simNWT.mat'];
    temp=load(simResFileName);simRes=temp.simRes;
    
    %Find time and location of steepest wave
    switch simRes.nwtSpec.solver.type
        case 'melhos'
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
                subplot(length(tWind),1,iSP)
                [~,indt]=min(abs(simRes.t-(tslMax+tWind(iSP))));
%                 plot(simRes.alpha+simRes.x(indt,:),simRes.z(indt,:),cCol{iSim})
                plot((simRes.alpha+simRes.x(indt,:))*km,simRes.z(indt,:)*km,[cCol{iSim} '-o'],'Markersize',1,'MarkerFaceColor',cCol{iSim})
            end
            
        case 'hos'
            [xMat,tMat]=meshgrid(simRes.x,simRes.t);
            dx=simRes.x(2)-simRes.x(1);
            slope=([simRes.eta(:,2:end) simRes.eta(:,1)]-[simRes.eta(:,end) simRes.eta(:,1:end-1)])./(2*dx);
            [slmax,indslMax]=max(-slope(:));
            if iSim==1
                indtslMax=find(simRes.t==tMat(indslMax));
                indxslMax=find(simRes.x==xMat(indslMax));
                xslMax=simRes.x(indxslMax);
                tslMax=simRes.t(indtslMax);
                xmin=xslMax+xWind(1);
                xmax=xslMax+xWind(2);
            end
            for iSP=1:length(tWind)
                subplot(length(tWind),1,iSP)
                [~,indt]=min(abs(simRes.t-(tslMax+tWind(iSP))));
                plot(simRes.x*km,simRes.eta(indt,:)*km,cCol{iSim},'linewidth',1)
            end
    end
    
end

%Finalize plot
for iSP=1:length(tWind)
    subplot(length(tWind),1,iSP)
    grid on
    axis equal
    title(['$t/T_m=' num2str((tslMax+tWind(iSP))/Tm,'%6.2f') '$' ],'interpreter','latex')
    ylabel('$k_m \hspace{0.05cm} z$','interpreter','latex')
    if iSP==1
        legend('MELHOS','HOS','Location','NorthWest')       
    end
    if iSP==length(tWind)
        xlabel('$k_m \hspace{0.05cm} x$','interpreter','latex');
    else
        set(gca,'XTickLabel','');
    end
    set(gca,'FontSize',9)
    set(gca,'XLim',[km*xmin km*xmax])
    ylim([-0.3 0.5])
end
saveas(f,'.\..\rapPlots\CompHOSMELHOS_O3','png')
set(gcf,'renderer','painters')
print(f,'.\..\..\..\Paper\ComputationalPhysics\melhos\CompHOSMELHOS_O3','-depsc')
