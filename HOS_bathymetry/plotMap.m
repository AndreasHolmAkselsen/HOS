function [hf_map,hf_mapZoom] = plotMap(fIP,xxIP,xxIP_near,xxLR,yyUpper,x,h0,zzS0,zzRoots)

hf_map = figure('color','w','position',[436 63 637 600]);

% plot the z-plane
haz = subplot(211); hold on;
%     haz = axes; hold on;
title('z-plane');xlabel('x');ylabel('i y');box off
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])

%     xxPlot = linspace(xxLR(1),xxLR(2),200);
%     xxPlot = xxIP;
%     xxPlot = xxIP(xxIP>=xxLR(1)&xxIP<=xxLR(2));
nxxPlot = 1000;
xxPlot = xxIP(1:round(length(xxIP)/nxxPlot):end);
%     zPhi = fzIP0(repmat(linspace(xxLR(1),xxLR(2),20)',1,200),repmat(linspace(-pi,yyUpperPlot,200),20,1)).';
%     zPsi = fzIP0(repmat(linspace(xxLR(1),xxLR(2),200)',1,10),repmat(linspace(-pi,yyUpperPlot,10),200,1)).';
zPhi = fIP(linspace(xxLR(1),xxLR(2),20)+1i*linspace(-pi,yyUpper,200).');
zPsi = fIP(xxPlot+1i*linspace(-pi,yyUpper,10).');
plot(zPhi,'r','linewidth',1); hold on; plot(zPsi.' ,'b')
minIz = min(imag(zPsi(1,:)));
patch(real(zPsi(1,[1,1:end,end])),[1.1*minIz,imag(zPsi(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');%,'FaceAlpha',.5
xlabel('x');ylabel('i y');
plot(x,h0,'k','linewidth',1.5);

% plot the zz-plane
hazz = subplot(212); hold on
title('\zeta-plane'); xlabel('\xi');ylabel('i \sigma'); box off
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
zzPlot = linspace(xxLR(1),xxLR(2),100) + 1i*linspace(-pi,yyUpper,100)';
zPlot = fIP(zzPlot);
contour(real(zzPlot),imag(zzPlot),real(zPlot),20,'r','linewidth',1);
contour(real(zzPlot),imag(zzPlot),imag(zPlot),10,'b','linewidth',1);
plot(zzS0,'k','linewidth',1.5);


% repeat for a 1-to-1 plot
%     xxLR = [-5,10];
xxLRnear = xxIP_near([1,end]);
hf_mapZoom = figure('color','w','position',[436 63 637 600]);
% plot the z-plane
haz = subplot(211); hold on;
%     haz = axes; hold on;
title('z-plane');xlabel('x');ylabel('i y');box off
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
zPhi = fIP(linspace(xxLRnear(1),xxLRnear(2),20)+1i*linspace(-pi,yyUpper,nxxPlot).');
zPsi = fIP(linspace(xxLRnear(1),xxLRnear(2),nxxPlot)+1i*linspace(-pi,yyUpper,10).');
plot(zPhi,'r','linewidth',1); hold on; plot(zPsi.' ,'b','linewidth',1)
minIz = min(imag(zPsi(1,:)));
patch(real(zPsi(1,[1,1:end,end])),[1.1*minIz,imag(zPsi(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');%,'FaceAlpha',.5
xlabel('x');ylabel('i y');
plot(x,h0,'k','linewidth',1.5);
plot(fIP(zzRoots+.025i),'ro'); % mark singularities
axis(haz,'equal','tight');xlim(haz,real(fIP(xxLRnear)));


% plot the zz-plane
hazz = subplot(212); hold on
title('\zeta-plane'); xlabel('\xi');ylabel('i \sigma'); box off
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
zzPlot = linspace(xxLRnear(1),xxLRnear(2),100) + 1i*linspace(-pi,yyUpper,100)';
zPlot = fIP(zzPlot);
contour(real(zzPlot),imag(zzPlot),real(zPlot),20,'r','linewidth',1);
contour(real(zzPlot),imag(zzPlot),imag(zPlot),10,'b','linewidth',1);
plot(zzS0,'k','linewidth',1.5);
plot(zzRoots,'ro'); % mark singularities
axis(hazz,'equal');xlim(hazz,xxLRnear);