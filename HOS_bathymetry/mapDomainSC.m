function [map,fIP,varphiS0,eta0_xiReg] = mapDomainSC(map,x,h0,phiS0,PLOT)

% Create inverse function through scattered interpolation:
for iNan = find(isnan(map.xx_b))
    map.xx_b(iNan) = map.xx_b(iNan-1) - pi./map.theta(iNan).*log(map.H(iNan+1)./map.H(iNan));
end
crudeScale = 1.5*pi/min(map.H);
yyUpper = max([2*max(h0)*crudeScale,.1*map.H]);
% xxIP = linspace(xLR(1),xLR(2),nArrX)*crudeScale;

xLR = [x(1),-x(1)];
xxIP_far = linspace(xLR(1),xLR(2),map.nArrX_far)*crudeScale;
map.zzRoots = [map.xx_b-1i*pi;map.xx_b+pi./map.theta.*log(map.H(2:end)./map.H(1:end-1))-1i*pi].';

xxRoots = real(map.zzRoots);
% xxIP_near = linspace( min(real(xxRoots(1,:)))-20*abs(diff(xxRoots(1,:))), max(real(xxRoots(end,:)))+5*abs(diff(xxRoots(end,:))), map.nArrX_near);
dxxRoot = max(xxRoots(:))-min(xxRoots(:));

if dxxRoot>0
    xxIP_near = linspace( min(real(xxRoots(1,:)))-5*dxxRoot, max(real(xxRoots(end,:)))+5*dxxRoot, map.nArrX_near);
else
    %         xxIP = xxIP_far;
    xxIP_near = 0; map.zzRoots( isnan(map.zzRoots))=0;
end
xxIP = [xxIP_far(xxIP_far<xxIP_near(1)), xxIP_near, xxIP_far(xxIP_far>xxIP_near(end))];

% yyIP = linspace(-pi,1.4*max(h0)*crudeScale,map.nArrYDown)';
yyIP = linspace(-pi,0,map.nArrYDown)'; % ensure that we capture the line yy=0
dy = yyIP(2)-yyIP(1);
yyIP = [yyIP; (dy:dy:yyUpper)' ];
assert(yyIP(map.nArrYDown)==0)

H = shiftdim(map.H,-1); theta = shiftdim(map.theta,-1); xx_b = shiftdim(map.xx_b,-1);
[zzIP,dfIP,zIP] = fmap_SCnum(xxIP,yyIP,H,theta,xx_b);
assert(real(zIP(map.nArrYDown,1))<x(1)&&real(zIP(map.nArrYDown,end))>x(end),'Physical domain [%.3g,%.3g] out of range of interpolation range [%.3g,%.3g]. Extend interpolation range.',x(1),x(end),real(zIP(map.nArrYDown,1)),real(zIP(map.nArrYDown,end)))

% plot interpolation basis to inspect resolution:
if PLOT
    figure('color','w');hold on
    plot(zIP(:,1:round(end/30):end),'r','linewidth',1);  plot(zIP(1:round(end/10):end,:).' ,'b')
    minIz = min(imag(zIP(1,:))); patch(real(zIP(1,[1,1:end,end])),[1.1*minIz,imag(zIP(1,:)),1.1*minIz],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
end

if ~all(abs(imag(zIP(map.nArrYDown,:)))<1e-3*max(map.H))
    warning('Map appears not to have a flat surface at yy=0. max|y(yy=0)| = %g',max(abs(imag(zIP(map.nArrYDown,:)))))
end
% assert(all(abs(imag(zIP(map.nArrYDown,:)))<1e-3*max(map.H)),'Map appears not to have a flat surface at yy=0')
xxLR = interp1(real(zIP(map.nArrYDown,:)),xxIP,xLR);

% iTrim = xxIP>=xxLR(1)&xxIP<=xxLR(2);
iTrim = (find(xxIP<xxLR(1),1,'last')-1):(find(xxIP>xxLR(2),1,'first')+1);
zzIP = zzIP(:,iTrim); zIP = zIP(:,iTrim); dfIP = dfIP(:,iTrim); %xxIP = xxIP(:,iTrim);
assert(all(real(zzIP(:,1))<xxLR(1)) && all(real(zzIP(:,end))>xxLR(2)))


fzIP0 = griddedInterpolant(real(zzIP).',imag(zzIP)',zIP.','linear','none');
fy0 = griddedInterpolant(real(zzIP).',imag(zzIP).',imag(zIP).','linear','none');
%     fJInv0 = griddedInterpolant(real(zzIP).',imag(zzIP).', abs(dfIP).'.^(-2) ,'linear','none');
f_zz0 = griddedInterpolant(real(zzIP).',imag(zzIP).', dfIP.','linear','none');

complex2grid = @(zz,ff) ff( (real(zz)+0*zz).', (imag(zz)+0*zz).').';
fIP = @(zz) complex2grid(zz,fzIP0);
map.fy = @(zz) complex2grid(zz,fy0);
%     map.fJInv = @(zz) complex2grid(zz,fJInv0);
map.f_zz =  @(zz) complex2grid(zz,f_zz0);
map.xxLR = xxLR;

finvIp = scatteredInterpolant(  real(zIP(:)), imag(zIP(:)) , zzIP(:),'linear','none');
zzS0 = finvIp(x, h0 );
% interpolate onto regulart xi grid.
map.xi = linspace(real(zzS0(1)),real(zzS0(end)),length(zzS0))';
map.zzDepth = pi;

eta0_xiReg = interp1(real(zzS0),imag(zzS0),map.xi);
xS_xiReg = real(fIP(map.xi+1i*eta0_xiReg));
if x(end) < xLR(2)-1e-12  %strcmp(boundaryType,'open')
    varphiS0 = interp1( [x-Lx;x;x+Lx],[phiS0;phiS0;phiS0],xS_xiReg );
else
    assert(x(1)<=xS_xiReg(1)&&x(end)>=xS_xiReg(end))
    varphiS0 = interp1( x,phiS0,xS_xiReg);
end
%     figure, contourf(real(zIP),imag(zIP),map.fJInv(zIP));colorbar
%     figure, subplot(121); plot(xS_xiReg,map.fy(map.xi+1i*eta0_xiReg),x,h0,'--');
%     subplot(122);plot(xS_xiReg,varphiS0,x,phiS0,'--');





