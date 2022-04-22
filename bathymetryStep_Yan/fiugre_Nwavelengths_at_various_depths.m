g=9.81;
L = 50;

h = 6:-.1:.25;

Nlambda = (2:5)'; % number of wave lengths on shoal

lambda = L./Nlambda;
k = 2*pi./lambda;
w = sqrt(k.*g.*tanh(k.*h));
T = 2*pi./w;
T_deep = T;
T(k.*h>pi/2)=nan;
T_deep(k.*h<pi/2) = nan;


% dw: deep water limit
hk_dw = [pi/2;pi];
h_dw = .1:.01:6;
k_dw = hk_dw./h_dw;
w_dw = sqrt(k_dw.*g.*tanh(hk_dw));
T_dw = 2*pi./w_dw;

figure('color','w');
plot(h,T,'-',h,T_deep,':','linewidth',1); grid on; hold on
hp = plot(h_dw,T_dw,'k--','linewidth',.5);
% hp(2).LineStyle = '-.';
text(4,2,'kh > \pi','fontsize',14)
text(2,2,'kh > \pi/2','fontsize',14)

xlabel('Water depth [m]');
ylabel('Period [s]')
ylim([0,14])
legend(Nlambda+" wavelengths")


return
savefig(['./memo/figures/fiugre_Nwavelengths_at_various_depths.fig']);
export_fig(['./memo/figures/fiugre_Nwavelengths_at_various_depths'],'-pdf');