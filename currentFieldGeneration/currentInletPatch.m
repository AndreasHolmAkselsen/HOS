clear
close all 

DO_EXPORT = true;
exportPrefix = 'basin_';
exportPath = './HOS_SFo_curr/basinCurrfigures/';

h = 15;
d = 4.5;
L = 50;

% nk = 512;
% % nk = 256;
% nx = 200;
% nz = 100;

nk = 128;
nx = 100;
nz = 50;


U0 = .17;
% VI = .5;
dk = pi/h;
k = shiftdim(1:nk,-1)*dk;

dz = h/nz;
z = (-nz+1:1:0)'*dz;
x = linspace(0,L,nx);

inletType = 'uniform';

Lambda = .5*k*h + .25*sin(2*k*h);
% LambdaS = .5*k*h - .25*sin(2*k*h);
switch inletType
    case 'uniform'
        UI =  U0/(1-d/h);
        Gamma = - UI.*sin(k*d)./k;    
%         GammaS = VI./k.*(cos(k*h)-cos(k*d));
    case 'linear'
        UI = 2*U0./(1-d/h);
        Gamma =  UI.*(cos(k*d)-cos(k*h)-(h-d)*k.*sin(k*d))./(k.^2*(h-d));
    case 'flipLinear'
        UI = 2*U0./(1-d/h);
        Gamma = UI.*(cos(k*h)-cos(k*d))./(k.^2.*(h-d));
end
hphi = -Gamma./Lambda;
% hphi = -Gamma./Lambda -GammaS./LambdaS;





% Current
zeta_j = 8-2.1i ;% object centre
U_j    =  -.55i ;% object strength-- +1:source, -1:sink, 1i:counter-clockwise vortex, -1i ...
patchWidth = 15;
patchHeight = 3.5;
nVortexes = 2*2*[5,20];

% patchWidth = 0;
% patchHeight = 0;
% nVortexes = [1,1];


zzPathch = linspace(0,patchWidth,nVortexes(2))- patchWidth/2 + 1i*(linspace(0,patchHeight,nVortexes(1)).'-patchHeight/2);
% zzPathch( imag(zzPathch)<(real(zzPathch)/patchWidth-.5)*patchHeight+patchHeight/2 )=[];
zzPathch = zeta_j + zzPathch(:).';
U_j = shiftdim(U_j,-1); zzPathch = shiftdim(zzPathch,-1);% ID_j = shiftdim(ID_j,-1);
% zeta_pathch = zeta_pathch + L*shiftdim(-nMirror:nMirror-1,-2);
A_j = .5*U_j.*abs(imag(zeta_j))/numel(zzPathch);
% f00  = @(zz) sum(A_j.*log(zz-zzPathch) + conj(A_j.*log(conj(zz)-zzPathch)),3:4);
% df00 = @(zz) sum(A_j./(zz-zzPathch) + conj(A_j./(conj(zz)-zzPathch)),3:4) ;
f00  = @(zz) sum(A_j.*log(zz-zzPathch) - A_j.*log(zz-conj(zzPathch)),3:4);
df00 = @(zz) sum(A_j./(zz-zzPathch) - A_j./(zz-conj(zzPathch)),3:4) ;
f0 = @(zz) f00(zz)   + f00(-zz);
df0 = @(zz) df00(zz) - df00(-zz);

% repeat entire domain once above and twice below
f = @(zz) f0(zz) + f0(zz+2i*h) + f0(zz-2i*h)  + f0(zz+4i*h);
df = @(zz) df0(zz) + df0(zz+2i*h) + df0(zz-2i*h)  + df0(zz+4i*h);



% g0  = @(zz) sum(A_j.*log(zz-zzPathch),3:4);
% dg0 = @(zz) sum(A_j./(zz-zzPathch),3:4) ;
% g = @(zz) g0(zz)+g0(-zz) + conj(g0(conj(zz))+g0(-conj(zz)));
% dg = @(zz) dg0(zz)-dg0(-zz) + conj(dg0(conj(zz))-dg0(-conj(zz)));


hfu = figure('color','w','position',[1640 694 560 284]);
u_x0 =  sum( -k.*hphi.*cos(k.*z),3 ) + U0 + real(df( 1i*z ));
w_x0 = sum( -k.*hphi.*sin(k.*z),3 ) - imag(df( 1i*z )) ;
subplot(1,2,1);plot(u_x0,z,'k',w_x0,z,'r');legend('u','w');title('Inlet velocity');ylabel('z');grid on
u_z0 =  sum( -k.*hphi.*exp(-k.*x),3 ) + U0 + real(df(x));
subplot(1,2,2);plot(x,u_z0,'k');title('Surface velocity');xlabel('x');grid on

hfc = figure('color','w','position',[242 406 1958 572]); hold on;
psi =  -sum( hphi.*sin(k.*z).*exp(-k.*x),3 ) + z.*U0 + imag(f(x+1i*z));
u =  -sum( k.*hphi.*cos(k.*z).*exp(-k.*x),3 ) + U0 + real(df(x+1i*z));
% [~,hc] = contourf(x,z,u,-.1:.01:.15,'linestyle','none');
[~,hc_u] = contourf(x,z,u,25,'linestyle','none');
phiLims = [min(psi(:)),max(psi(:))];
psi_middle = -phiLims(2);
[~,hc_psi] = contour(x,z,psi,[linspace(phiLims(1),psi_middle,15),linspace(psi_middle,phiLims(2),10)],'k');
caxis(hc_u.LevelList([1,end]))
colorbar;
zq = linspace(-h,0,10); 
u0q =  sum( -k.*hphi.*cos(k.*zq),3 ) + U0 + real(df( 1i*zq ));
w0q =  sum( -k.*hphi.*sin(k.*zq),3 ) - imag(df( 1i*zq )) ;
quiver(0*zq,zq,u0q,w0q,'r')
axis equal 
xlabel('x [m]','fontsize',12)
ylabel('z [m]','fontsize',12)
ylim([-h,0])
plot(zzPathch(:),'.k','MarkerSize',3)


if DO_EXPORT
    if ~isfolder(exportPath), mkdir(exportPath); end
    exportName = sprintf('%s%s_nVort%dx%d_dimVort%gx%g_Uj%g',exportPrefix,inletType,size(nVortexes),patchWidth,patchHeight,abs(U_j)); exportName(exportName=='.')='p';
    savefig(hfu,[exportPath,'/u_',exportName]);
    export_fig(hfu,[exportPath,'/u_',exportName],'-pdf','-png');    
    savefig(hfc,[exportPath,'/',exportName]);
    export_fig(hfc,[exportPath,'/',exportName],'-pdf','-png');
end