clear
h = 15;
T = .1:.01:2;
UR = 0.17;
UL = -0.05;

addpath C:\gits\timsas2\matlabLibs\
kR = level1.wave.computeWaveNumber(2*pi./T,h,UR,0);
kL = level1.wave.computeWaveNumber(2*pi./T,h,UL,0);

lamR = 2*pi./kR;
lamL = 2*pi./kL;

iRm = imag(lamL)~=0;
T(iRm)=[];lamR(iRm)=[];lamL(iRm)=[];

figure('color','w','position',[1207 612 481 367]); 
hold on; grid on; box on;
xlabel('Period [s]'); ylabel('Wavelength [m]');
plot(T,lamL,T,lamR,T,lamL./lamR,'linewidth',1.5);
legend({'\lambda_{min}','\lambda_\infty','\lambda_{min}/\lambda_\infty'},'location','northwest','fontsize',11)
axis tight
ylim([0,3])

return
export_fig('./doc/figures/lambda_envelope','-pdf')