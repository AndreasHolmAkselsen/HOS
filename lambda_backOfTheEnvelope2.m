clear
h = 10;
VelMultiplyer = (.5:.5:2)';
UR0 = 0.17;
UL0 = -0.05;
g = 9.81;
T = linspace(0,5,100);

addpath C:\gits\timsas2\matlabLibs\
% for i=1:length(h)
% kR(i,:) = level1.wave.computeWaveNumber(2*pi./T,h(i),UR,0);
% kL(i,:) = level1.wave.computeWaveNumber(2*pi./T,h(i),UL,0);
% end
w=2*pi./T;
UL = UL0*VelMultiplyer;UR = UR0*VelMultiplyer;
for i=1:length(VelMultiplyer)
    kR(i,:) = level1.wave.computeWaveNumber(w,h,UR(i),0);
    kL(i,:) = level1.wave.computeWaveNumber(w,h,UL(i),0);
end
kL(imag(kL)~=0) = nan;
% irm = imag(kL)~=0;
% kR(:,irm)=[];kL(:,irm)=[];T(:,irm)=[];

% cR = sqrt(g./kR.*tanh(kR*h));
% T = (1+UR./cR).*tanh(kR*h)./tanh(kL*h);
% cL__cR = .5./T.*(1+sqrt(1+4*T.*UL./cR));
% cL = cL__cR.*cR;
% cL_test =  sqrt(g./kL.*tanh(kL*h));
% figure, plot(T,cL,'-',T,cL_test,'--','linewidth',1.5)

% HbL = .*.88./kL.*tanh(.89*kL*h); % Miche's formula

% AL = HbL/2;
% epsL = AL.*kL;

sigL = w-UL.*kL; sigR = w-UR.*kR;
cgL = .5*sigL./kL.*(1+2*kL.*h./sinh(2*kL.*h));
cgR = .5*sigR./kR.*(1+2*kR.*h./sinh(2*kR.*h));
AR__AL = sqrt( (UL+cgL)./sigL .* sigR./(UR+cgR) );

epsR__epsL = AR__AL.*kR./kL;

figure('color','w','position',[1207 612 481 367]);
hold on; grid on; box on;
xlabel('Period [s]'); ylabel('$(ka)_\infty/(ka)_\mathrm{min}$','fontsize',14,'interpreter','latex');
plot(T,epsR__epsL,'linewidth',1.5);
% title("$U_\mathrm{min}^0="+UL0+"\,$m/s$, U_\infty^0="+UR0+"\,$m/s",'interpreter','latex');
legend("$(U_\mathrm{min},U_\infty) = "+VelMultiplyer+"\times(U_\mathrm{min}^0,U_\infty^0)$",'location','southeast','interpreter','latex','fontsize',11)
% axis tight

return
export_fig('./doc/figures/waveActionSteepness','-pdf','-png')