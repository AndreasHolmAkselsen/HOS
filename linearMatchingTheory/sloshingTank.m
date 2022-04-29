clear
close all
g = 9.81;
addpath c:/gits/wavemaker

L = 5;
T = 1;
w = 2*pi/T;
h = 1;
nEv = 200;
t = 100*rand;
U00 = .04;
x = linspace(0,L,100);
z = linspace(-h,0,100).';

sourceType = 'grid';  d = .135*h;
sourceType = 'heaviside';  d = .4*h;
sourceType = 'sin';  d = .2*h;
sourceType = 'U0*z';  d = .4*h;
%% code

bad=@(x) any(isnan(x)|isinf(x),'all');

k = findWaveNumbers(w,h,0,nEv);
k = shiftdim([k;-k],-2);

U0 = U00*2*pi/k(1)/T;
kh = k*h;
switch sourceType
    case 'U0*z'
        I = U0./((h-d)*k.^2).*((cosh(k*(h-d))-1)./cosh(kh));
    case 'sin'
        I = U0*pi*d./(pi+(k*d).^2).*(cos(pi*h/d).*sech(kh)-1);
    case 'heaviside'
        I = U0./k.*sinh(k*(h-d)).*sech(kh);
    case 'grid'
        N = ceil(.5*h/d);
        j = 1:N; % z_j = -min(j*d,h);
        I = U0./k.*sech(kh).*sum( sinh(k.*(h-(2*j-2)*d)) - sinh(k.*(h-(2*j-1)*d)) ,2);
end
kLambda = .5*(kh.*sech(kh).^2+tanh(kh));

frontTerm =  .5*exp(-1i.*k.*L)./sin(k*L);
frontTerm(imag(k)*L>100) = -1i; % limit
hphi = frontTerm ./kLambda .* I;


exp_ikx =  exp(1i*k.*x);
exp_ikx(:,:,hphi==0) = 0;
phi0 = sum(  hphi.*cosh(k.*(h+z))./cosh(kh).*exp_ikx,3);
phi = phi0.*cos(w*t);
eta = w*phi0(1,:)/g.*sin(w*t);

phix_z0 = sum( 1i*k.* hphi.*cosh(k.*(h+z))./cosh(kh),3);


hf = figure('color','w','position',[1640 558 1129 420]);
subplot(1,2,1);
contourf(x,z,real(phi));
hold on
plot(x,eta,'k','linewidth',1)

ha=subplot(1,2,2);
plot(phix_z0,z,'k','linewidth',1);
ylabel('z'); grid on;
ha.XAxisLocation='origin';
xlabel('U(z)');

