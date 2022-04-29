clear all
close all
PLOT_OTHER_PLANES = 0;

wbl = .6;
wbOveWater = .1;
thetaMax = 20*pi/180;
L = .4;

limZeta = [0.001,.7*wbl,-.5*wbl,.5*wbl];
% limZeta = [0.001,2,-1,1];

cLines = 10;
% cLines = (0:.1:1)*.5*wbl;
nPoints = 100;%500;

dt = 0.1;
tMax = 10;
T = 1;
nTerms = 50;


% d = 1; % unit domain 
% d = wbl + wbOveWater; % scale two domains


zeta = linspace(limZeta(1),limZeta(2),nPoints) + 1i*linspace(limZeta(3),limZeta(4),nPoints)';

hf = figure('color','w');
hold on;
t = 0;
k = 1;
while t<tMax
    theta = thetaMax*cos( 2*pi/T*t );
    thp = theta/pi;

%     H = hypergeometric2f1(1,1.5,1.5-thp,-(zeta/d).^2,100);
%     fz = -tan(theta).*(d+gamma(-thp)./gamma(.5-thp).*sqrt(pi)/(pi-2*theta).*zeta.^(1-2*thp).*(zeta-1i*d).^thp.*(zeta+1i*d).^thp.*(1+(zeta/d).^2).*H);

    c = -(wbl+wbOveWater)*tan(theta) + 1i*wbOveWater;
%     k = pi^1.5/(pi-2*theta).*(wbl-1i*c).*exp(-1i*theta)./(gamma(1+thp)*gamma(.5-thp));
%     d = L*real(k)*(L/wbl)^(-2*thp)*hypergeometric2f1(-thp,.5-thp,1.5-thp,-(L/wbl).^2,nTerms)/(L+wbl*tan(theta));   
%     d = k*(1-2*thp);

    d = sqrt(pi).*(wbl-1i*c).*exp(-1i*theta)./(gamma(1+thp)*gamma(.5-thp));
    fz = c + k*(zeta/d).^(1-2*thp).*hypergeometric2f1(-thp,.5-thp,1.5-thp,-(zeta/d).^2,nTerms);
    
    cla
    contour(real(fz),imag(fz),real(zeta),cLines,'r','linewidth',1);
    contour(real(fz),imag(fz),imag(zeta),cLines,'b','linewidth',1);
    t = t+dt;
    drawnow
end

% dat = load('FF_r.csv')+1i* load('FF_i.csv')


