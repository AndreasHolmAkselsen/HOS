clear
close all 


DO_PLOT = true;
DO_EXPORT_FIG = true;
DO_EXPORT_MAT = true;
exportPath = './currendDatabase/';
exportPrefix = '';
includeZPositive = true;


% % read from existing...
% temp=readstruct(['.\basinSpec\currentTestSmall.specNWT']);
% nwt=temp.nwt;
% or specify here
nwt.basin.length =  84.8694;% 200;
nwt.basin.depth = 15;
nwt.current.wallDepth = 0;%4.5;
nwt.current.nk = 512;
nwt.current.U0 = .17;
nwt.current.inletType = 'uniform';
nwt.vortexPatch.nx = 1024;
nwt.vortexPatch.zRange0 = 2;
% nwt.vortexPatch.patchCentre = 8-2.1i;
% nwt.vortexPatch.patchCentre = .5*nwt.basin.length -1.8i;
nwt.vortexPatch.patchCentre = 50 -1.8i;
nwt.vortexPatch.patchStrength = -.55i;
nwt.vortexPatch.patchSize = [15,3.5];

% nwt.vortexPatch.patchStrength = 0;
% nwt.vortexPatch.patchSize = [0,0];

nwt.vortexPatch.mirror = false;

% % non-staggered field:
% nwt.vortexPatch.nx = 4*1024;
% nwt.vortexPatch.patchCentre = 8-2.5i;
% nwt.vortexPatch.nSize = [10,40]; % if present, then the staggered grid is overwritten

% % for full basin picture:
% nwt.vortexPatch.zRange0 = 15;
% nwt.vortexPatch.nx = 256;
% includeZPositive = false;
% nwt.vortexPatch.mirror = true;

%% code

nwt = current.setCurrentParam(nwt);
nx = nwt.vortexPatch.nx;
patchCentre = nwt.vortexPatch.patchCentre;
patchStrength = nwt.vortexPatch.patchStrength;
patchSize = nwt.vortexPatch.patchSize;

[dx,dz] = deal(nwt.basin.length/nx); % uniform grid!
nz = round( nwt.vortexPatch.zRange0/dz );
zrange = nz*dz; % adjust
z = ( -nz+1:includeZPositive*(nz-1) )'*dz;
x = (0:nx-1)*dx;

if isfield(nwt.vortexPatch,'nSize')
    % old method
    zzPatch = nwt.vortexPatch.patchCentre + linspace(0,nwt.vortexPatch.patchSize(1),nwt.vortexPatch.nSize(2))-nwt.vortexPatch.patchSize(1)/2 + 1i*(linspace(0,nwt.vortexPatch.patchSize(2),nwt.vortexPatch.nSize(1)).'-nwt.vortexPatch.patchSize(2)/2);
    exportPrefix = sprintf('%snonStag_nSize%dx%d',exportPrefix,nwt.vortexPatch.nSize);
else
    % new method: fit to structured grid.
    jc = round( (real(patchCentre)+dx/2-x(1))/dx );
    ic = round( (imag(patchCentre)+dz/2-z(1))/dx );
    xVortex = x(1)+(( -round(patchSize(1)/2/dx):round(patchSize(1)/2/dx)  )+jc-.5)*dx;
    zVortex = z(1)+(( -round(patchSize(2)/2/dz):round(patchSize(2)/2/dz)  )+ic-.5).'*dz;
    % figure;
    % plot(x,z,'k.',xVortex,zVortex,'rx','markersize',4);hold on
    % plot(patchCentre,'b+','markersize',10);
    % plot(patchCentre+.5*[-patcSize(1);patcSize(1);-patcSize(2);patcSize(2)],'g^')
    % axis equal
    zzPatch = xVortex + 1i*zVortex;
end

% zzPatch( imag(zzPatch)<(real(zzPatch)/width-.5)*height+height/2 )=[];
nwt.vortexPatch.zzPatch = shiftdim(zzPatch(:),-2);
nwt.vortexPatch.A_j = .5*patchStrength.*abs(imag(patchCentre))/numel(nwt.vortexPatch.zzPatch);


exportName = sprintf('%s%s_Lh%gx%g_nx%dnk%d_patchPos%gx%gx%gx%g_Uj%g_posZ%d',exportPrefix,nwt.current.inletType,nwt.basin.length,nwt.vortexPatch.zRange0,nx,nwt.current.nk,real(nwt.vortexPatch.patchCentre),imag(nwt.vortexPatch.patchCentre),patchSize,abs(patchStrength),includeZPositive); 
exportName(exportName=='.')='p';exportName(exportName=='-')='m';
if DO_PLOT
    [u,w,phi,psi] =  current.velocityFieldVortexPatch(x,z,nwt);    
    
%     u_x0 = u(:,x==0); w_x0 = w(:,x==0); z_z0=z;
    z_z0 = 0:-dz:-nwt.basin.depth;
    [u_x0,w_x0] = current.velocityFieldVortexPatch(0,z_z0,nwt);
    %     u_z0 = current.velocityFieldVortexPatch(x,0,nwt);
    
    hfu = figure('color','w','position',[2391 491 1154 347]); % [1640 694 560 284]
    subplot(1,3,1);plot(u_x0,z_z0,'k',w_x0,z_z0,'r');legend('u','w');title('Inlet velocity');ylabel('z');grid on
    subplot(1,3,2);plot(x,u(z==0,:),'k');title('Horizontal surface velocity');
    xlabel('x');ylabel('u'); grid on
    ii = find(abs(z)<.3);ii = ii(1:2:end);
    subplot(1,3,3);plot(x,w(ii,:));legend("z = "+round(z(ii),2)'+" m")
    title('Vertival velocity near surface')
    xlabel('x'); ylabel('w'); grid on;
    
    
    
    hfc = figure('color','w','position',[242 406 1958 572]); hold on;
    [~,hc_u] = contourf(x,z,u,25,'linestyle','none');
    phiLims = [min(psi(z+0*x<0)),max(psi(z+0*x<0))];
%     psi_contours = [linspace(phiLims(1),-phiLims(2),15),linspace(-phiLims(2),phiLims(2),10)];
    psi_contours = [linspace(phiLims(1),0,15)];
    [~,hc_psi] = contour(x,z,psi,psi_contours,'k');
    try,caxis(hc_u.LevelList([1,end]));colorbar;end
    % zq = linspace(-nwt.basin.depth,0,10);
    zq = linspace(z(1),z(end),10);
    [u0q,w0q] =  current.velocityFieldVortexPatch(0,zq,nwt);
    quiver(0*zq,zq,u0q,w0q,'r')
    axis equal
    xlabel('x [m]','fontsize',12)
    ylabel('z [m]','fontsize',12)
    ylim(z([1,end]))
    
    plot(zzPatch(:),'.k','MarkerSize',3)
    % [xx,zz]=meshgrid(x,z);plot(xx(:),zz(:),'.k','MarkerSize',1)
    
    if DO_EXPORT_FIG
        if ~isfolder(exportPath), mkdir(exportPath); end
        savefig(hfu,[exportPath,'/u_',exportName]);
        export_fig(hfu,[exportPath,'/u_',exportName],'-pdf','-png');
        savefig(hfc,[exportPath,'/',exportName]);
        export_fig(hfc,[exportPath,'/',exportName],'-pdf','-png');
    end
end


if DO_EXPORT_MAT
    [xx,zz]=ndgrid(x,z);
    nwt.current.dfInterp = griddedInterpolant(xx,zz,(u-1i*w).');
    % % plot some random lines:
    % [xt,zt] = ndgrid(0:.001:5,z(1)*rand(10,1));
    % figure; plot(xt,-imag(nwt.current.dfInterp(xt,zt)),'r',xt,real(nwt.current.dfInterp(xt,zt)),'b')
    
    if ~isfolder(exportPath), mkdir(exportPath); end
%     matExportName = specFileName;
    save([exportPath,'/',exportName],'nwt')
    save([exportPath,'/workspace_',exportName])
    
    % save struct as specNWT file.
    nwtExp = nwt;
    nwtExp.vortexPatch = rmfield(nwtExp.vortexPatch,'zzPatch');
    nwtExp.current=rmfield(nwtExp.current,{'k','hphi','dfInterp'});
    fid = fopen([exportPath,'/',exportName,'.specNWT'],'w');
    writestruct(fid,'nwtExp');
    fclose(fid);
end


return
hfw = figure('color','w','position',[1640 694 560 284]);
ii = find(abs(z)<.3);ii = ii(1:2:end);
plot(x,w(ii,:));legend("z = "+round(z(ii),2)'+"m")
title('Vertival velocity near surface')
xlabel('x'); ylabel('w');
