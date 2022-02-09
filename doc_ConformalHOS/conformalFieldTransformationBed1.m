close all

% k = -1;
% a = .1;
% L = 2*pi;
% h = 1;

k0 = pi;
A = .075*exp(1i*30*pi/180 );
L = 2*pi/k0;
H = .75;

heta = shiftdim( .75*(A+A') ,-1);
k = shiftdim( [1,-1]*k0  ,-1);



% z|-> f(zeta)
f =@(zeta) zeta - 1i*sum( heta.*exp(1i*k.*(zeta+1i*H))./sinh(k*H) ,3);

xi = linspace(-L,L,100);
xS = real(f(xi));
eta = imag(f(xi));
% equivalent to 
% eta = sum( heta.*exp(1i*k.*xi) ,3);
% xS =  xi - 1i*sum( heta.*exp(1i*k.*xi)./tanh(k*H) ,3)  


hf = figure('color','w');
% clf
[Xi,Sig]= meshgrid(linspace(-L,L,10),linspace(-H,.25*H,100));
plot(f(Xi+1i*Sig),'r'); hold on;

[Xi,Sig]= meshgrid(linspace(-L,L,100),linspace(-H,.25*H,10));
plot(f(Xi+1i*Sig).','b')
axis equal
grid on
axis tight

plot(xS,eta,'k','linewidth',1.5)
plot(xS([1,end]),-H*[1,1],'k','linewidth',1.5)
xlabel('x')
ylabel('i y')


return
export_fig(hf,'conformalH','-pdf')


