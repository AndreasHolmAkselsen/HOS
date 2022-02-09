clear variables

outPath='.\..\simResults';

%Read simulation spec
% specFile='testMELHOS_BichromWave2.1';
specFile='testMELHOS_Gaussian';
temp=readstruct(['.\..\inputFiles\',specFile,'.specNWT']);
nwt=temp.nwt;

%Run hos-simulation
simRes=runNWT(nwt);
save([outPath '\' specFile,'_simNWT.mat'],'simRes','-MAT');

figure;imagesc(simRes.t,simRes.alpha,simRes.eta.')
figure;plot(simRes.t,simRes.eta(:,length(simRes.alpha)/2));grid on

%Plot eta(x) for largest wave time
[alphaMat,tMat]=meshgrid(simRes.alpha,simRes.t);
[zmax,indzMax]=max(simRes.z(:));
indtMax=find(simRes.t==tMat(indzMax));
figure
hold on
plot(simRes.alpha,simRes.eta(indtMax,:),'b')
plot(simRes.alpha+simRes.x(indtMax,:),simRes.z(indtMax,:),'r--','linewidth',2')
grid on

%Plot eta for largest front steepness time
dalpha=simRes.alpha(2)-simRes.alpha(1);
slope=([simRes.z(:,2:end) simRes.z(:,1)]-[simRes.z(:,end) simRes.z(:,1:end-1)])./...
    (2*dalpha+[simRes.x(:,2:end) simRes.x(:,1)]-[simRes.x(:,end) simRes.x(:,1:end-1)]);
[slmax,indslMax]=max(-slope(:));
indtslMax=find(simRes.t==tMat(indslMax));
indaslMax=find(simRes.alpha==alphaMat(indslMax));
xslMax=simRes.alpha(indaslMax)+simRes.x(indtslMax,indaslMax);
tslMax=simRes.t(indtslMax);

% figure;imagesc(simRes.t,simRes.alpha,slope.');colorbar

% figure
% hold on
% plot(simRes.alpha+simRes.x(indtslMax,:),simRes.z(indtslMax,:),'b');
% axis equal
% xlim([xslMax-simRes.nwtSpec.sim.Lx/20 xslMax+simRes.nwtSpec.sim.Lx/20])
% grid on
% 
% figure;
% hold on
% plot(simRes.alpha+simRes.x(indtslMax,:),slope(indtslMax,:),'b');
% grid on


%Compute velocity field
% xslMax=78;
% xslMax=43.95;
xOut=(xslMax-10):0.05:(xslMax+10);
zOut=-1:0.01:zmax;
tOut=(tslMax-10):0.04:(tslMax+10);
[u,w,tUsed]=computeNWTVelocity(simRes,xOut,zOut,tOut);
save([outPath '\' specFile,'_velNWT.mat'],'xOut','zOut','tUsed','u','w','-MAT');

%Make video of velocity field
iSkip=1;
figure('Position',[0 0 1500 350])
v=VideoWriter([outPath '\VelocityField_u_',specFile,'_t',...
    num2str(min(tOut)),'_',num2str(max(tOut)),'.avi']);
v.FrameRate=fix(1/(tOut(2)-tOut(1)));
open(v);
for indT=1:iSkip:length(tUsed)
    iT=find(simRes.t==tUsed(indT));
    set(gca,'nextplot','replacechildren');
    imagesc(xOut,zOut,u(:,:,indT).','AlphaData',~isnan(u(:,:,indT).'))
    hold on
    plot(simRes.alpha+simRes.x(iT,:),simRes.z(iT,:),'k')
%     plot(simRes.alpha,simRes.eta(iT,:),'k','linewidth',2)
    plot(simRes.alpha(1:4:end)+simRes.x(iT,1:4:end),simRes.z(iT,1:4:end),'o','MarkerFaceColor','b','Markersize',3)
    plot(simRes.alpha(1:8:end)+simRes.x(iT,1:8:end),simRes.z(iT,1:8:end),'o','MarkerFaceColor','k','Markersize',6)
    set(gca,'YDir','normal')
    caxis([min(u(:)) max(u(:))])
    colorbar
    grid on
    xlabel('x [m]')
    ylabel('\eta(x)')
    xlim([min(xOut) max(xOut)])
    axis equal
    title(['Horizontal velocity [m/s] - ',specFile,'t=' num2str(tUsed(indT),'%6.2f'),'s'],'interpreter','none')
    set(gca,'FontSize',12)
    colormap jet
    
    frame=getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);
