close all



if ~exist('siU','var')
    csvFile = 'C:\Users\andreasak\SINTEF\OSC - Current system CFD simulations folder - Dokumenter\BasinDesignSep2021\CFD_results\4Csaba\Case000\XYZ_Internal_Table_Y=25_6000.csv';
    tab = readtable(csvFile);
    assert( numel(unique(tab.Y_m_))==1, 'data not constant in y' );
    siU = scatteredInterpolant(tab.X_m_,tab.Z_m_-max(tab.Z_m_),tab.Velocity_i__m_s_,'linear','none');
    siW = scatteredInterpolant(tab.X_m_,tab.Z_m_-max(tab.Z_m_),tab.Velocity_k__m_s_,'linear','none');
end

nx = 200;
nz = 100;

x = linspace(0,50,nx);
z = linspace(0,-1,nz).';


% [xx,zz] = meshgrid(linspace(min(x),max(x),200),linspace(min(z),max(z),100));
[xx,zz] = meshgrid(x,z);
Usurf = siU(x,max(z)+0*x);
Ulim = [min(Usurf),1.2*max(Usurf)];
hf = figure('color','w','position',[2503 345 1236 576]);
caxis(Ulim);
subplot(2,1,1); contourf(xx,zz,siU(xx,zz),linspace(Ulim(1),Ulim(2),20));hold on
colorbar;title('U')
axis equal;
subplot(2,1,2); contourf(xx,zz,siW(xx,zz),20);
colorbar; title('W')
axis equal;

[phi1,phi2,psi1,psi2] = deal(zeros(size(xx)));

xh = .5*(x(1:end-1)+x(2:end));
zh = .5*(z(1:end-1)+z(2:end));
dz = diff(z);
dx = diff(x);
phi1(:,1) = [0; cumsum(siW(x(1)+0*zh, zh).*dz,1)]; 
phi1(:,2:nx) = phi1(:,1) + cumsum(siU(xh+0*z,z+0*xh).*dx,2); 

psi1(:,1) = [0; cumsum(siU(x(1)+0*zh,zh).*dz,1)]; 
psi1(:,2:nx) = psi1(:,1) - cumsum(siW(xh+0*z,z+0*xh).*dx,2); 

phi2(1,:) = [0, cumsum(siU(xh,z(1)+0*xh).*dx,2)];
phi2(2:nz,:) = phi2(1,:) + cumsum(siW(x+0*zh,zh+0*x).*dz,1);

psi2(1,:) = [0, -cumsum(siW(xh,z(1)+0*xh).*dx,2)];
psi2(2:nz,:) = phi2(1,:) + cumsum(siU(x+0*zh,zh+0*x).*dz,1);

omega_j = nan(size(xx));
xxc = xx(2:nz-1,2:nx-1);zzc = zz(2:nz-1,2:nx-1);
omega_j(2:nz-1,2:nx-1) =  (siU(xxc,zz(3:nz,2:nx-1)) - siU(xxc,zz(1:nz-2,2:nx-1)))./(z(3:nz)-z(1:nz-2))...
                        - (siW(xx(2:nz-1,3:nx),zzc) - siW(xx(2:nz-1,1:nx-2),zzc))./(x(3:nx)-x(1:nx-2));

hf2 = figure('color','w','position',[2503 49 1236 872]);
subplot(3,1,1);hold on;
[~,hppji1] = contourf(xx,zz,phi1,20);
contour(xx,zz,psi1,'k');
caxis([min(phi1(:)),max(phi1(:))]);colorbar
subplot(3,1,2);hold on
% contourf(xx,zz,phi2,hppji1.LevelList);
contourf(xx,zz,phi2,20);
contour(xx,zz,psi2,'k');
caxis([min(phi2(:)),max(phi2(:))]);colorbar
subplot(3,1,3);hold on
contourf(xx,zz,omega_j,linspace(-.2,.2,20));colorbar;

ws = load('C:\gits\SFo_gitHOS\hosm-nwt2d\basinSpec\workspace_flipLinear_Lh50x1_nx1024nk512_dimVort15x3p5_Uj0p55_posZ1.mat');
figure('color','w','position',[2503 345 1236 576])
subplot(2,1,1); contourf(ws.xx,ws.zz,ws.u.',20);hold on
colorbar;title('U')
axis equal;
ylim([-1,0])
subplot(2,1,2); contourf(ws.xx,ws.zz,ws.w.',20);
colorbar; title('W')
axis equal;
ylim([-1,0])