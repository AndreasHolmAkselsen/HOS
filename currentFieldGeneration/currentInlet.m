clear

h = 2.2;
d = .15*h;

nk = 512;
nx = 200;
nz = 200;

U0 = .5;
% VI = .5;
dk = pi/h;
k = shiftdim(1:nk,-1)*dk;

dz = h/nz;
z = (-nz+1:1:0)'*dz;
x = linspace(0,h,nx);

inletType = 'flipLinear';

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
        Gamma = UI.*(cos(k*d)-cos(k*h))./(k.^2.*(d-h));
end
hphi = -Gamma./Lambda;
% hphi = -Gamma./Lambda -GammaS./LambdaS;

phi0_x = sum( -k.*hphi.*cos(k.*z),3 ) + U0;
phi0_z = sum( -k.*hphi.*sin(k.*z),3 ) ;
figure
plot(phi0_x,z,'k',phi0_z,z,'r');

figure
phi =  sum( hphi.*cos(k.*z).*exp(-k.*x),3 ) +x.*U0;
contourf(x,z,phi)
hold on;
zq = linspace(-h,0,10); 
phi0_xq = sum( -k.*hphi.*cos(k.*zq),3 ) + U0;
phi0_zq = sum( -k.*hphi.*sin(k.*zq),3 ) ;
quiver(0*zq,zq,phi0_xq,phi0_zq,'r')