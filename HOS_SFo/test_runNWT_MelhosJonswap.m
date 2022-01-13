clear variables

outPath='.\..\simResults';

%Read simulation spec
specFile='testMELHOS_RegWaveM6_eps0.1';
temp=readstruct(['.\..\inputFiles\',specFile,'.specNWT']);
nwt=temp.nwt;

%Run hos-simulation
simRes=runNWT(nwt);
save([outPath '\' specFile,'_simNWT.mat'],'simRes','-MAT');

%Compute velocity field
xOut=10:0.05:30;
zOut=-2:0.01:0.60;
tOut=20:0.04:30;
[u,w,tUsed]=computeNWTVelocity(simRes,xOut,zOut,tOut);
save([outPath '\' specFile,'_velNWT.mat'],'xOut','zOut','tUsed','u','w','-MAT');
% 
% %Make video of velocity field
% iSkip=1;
% figure('Position',[0 0 1500 350])
% v=VideoWriter([outPath '\VelocityField_u_',specFile,'_t',...
%     num2str(min(tOut)),'_',num2str(max(tOut)),'.avi']);
% v.FrameRate=fix(1/(tOut(2)-tOut(1)));
% open(v);
% for indT=1:iSkip:length(tUsed)
%     iT=find(simRes.t==tUsed(indT));
%     set(gca,'nextplot','replacechildren');
%     imagesc(xOut,zOut,u(:,:,indT).','AlphaData',~isnan(u(:,:,indT).'))
%     hold on
%     plot(simRes.alpha,simRes.eta(iT,:),'k','linewidth',2)
%     plot(simRes.alpha(1:4:end)+simRes.x(iT,1:4:end),simRes.z(iT,1:4:end),'o','MarkerFaceColor','b','Markersize',4)
%     plot(simRes.alpha(1:8:end)+simRes.x(iT,1:8:end),simRes.z(iT,1:8:end),'o','MarkerFaceColor','k','Markersize',8)
%     set(gca,'YDir','normal')
%     caxis([min(u(:)) max(u(:))])
%     colorbar
%     grid on
%     xlabel('x [m]')
%     ylabel('\eta(x)')
%     xlim([min(xOut) max(xOut)])
%     axis equal
%     title(['Horizontal velocity [m/s] - ',specFile,'t=' num2str(tUsed(indT),'%6.2f'),'s'],'interpreter','none')
%     set(gca,'FontSize',12)
%     colormap jet
%     
%     frame=getframe(gcf);
%     writeVideo(v,frame);
%     
% end
% close(v);
