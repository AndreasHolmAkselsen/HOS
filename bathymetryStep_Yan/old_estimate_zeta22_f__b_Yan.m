clear
h_s = logspace(log10(2),log10(9.5),40)';
f = logspace(log10(.1),log10(.4),45);
% h_s = (2:.5:9.5)';
% f = (.1:.01:.36);
h_d = 10;
g = 9.81;


addpath .\step_functions_YanLi



% [k0hd,hshd] = meshgrid([.3,.5,1,2,4,6],[.2,.4,.6,.8,1]);
% T0Ref=0*k0hd;
% T0Ref(:,1) = [1.34, 1.205, 1.13, 1.05,1.025];
% T0Ref(:,2) = [1.32,1.2,  1.12,1.05,1.04];
% T0Ref(:,3) = [1.2, 1.18, 1.1, 1.04,1.03];
% T0Ref(:,4) = [1.11, .95, .95,  .98, .99];
% T0Ref(:,5) = [.99, .93,  .95, .98,.98];
% T0Ref(:,6) = [ .98, .97,  .98, .98,.99];
% 
% 
% T20Ref=0*k0hd;
% eps = 1e-6;
% T20Ref(:,1) = [23,23,22,21,1];
% T20Ref(:,2) = [22,20,20,8,1];
% T20Ref(:,3) = [20, 5, 2,1,eps];
% T20Ref(:,4) = [3,   1,  .5,.2,eps];
% T20Ref(:,5) = [.5, .1,  eps,eps,eps];
% T20Ref(:,6) = [eps, eps,  eps,eps,eps];
% 
% 
% R0Ref=0*k0hd;
% eps = 1e-6;
% R0Ref(:,1) = [.35  ,.210,.125,.055,eps];
% R0Ref(:,2) = [.33,.208,.12,.05,eps];
% R0Ref(:,3) = [.305, .195,.105,.045,eps];
% R0Ref(:,4) = [.215,.15,.06,.02,eps];
% R0Ref(:,5) = [.1,  .045,.02,eps, eps];
% R0Ref(:,6) = [.05, .01,eps,eps,eps];


w0 = 2*pi*f;
k0 = findWaveNumbers(w0,h_d,0,0);
% k0s = 0*w0.*h_s;
% for i = 1:length(h_s)
%     k0s(i,:) = findWaveNumbers(w0,h_s(i),0,0);
% end

fill = zeros(numel(h_s),numel(f));

% R0 = interp2(k0hd,hshd,R0Ref,  k0.*h_d+fill,  h_s./h_d+fill,'makima',nan);
% T0 = interp2(k0hd,hshd,T0Ref,  k0.*h_d+fill,  h_s./h_d+fill,'makima',nan);
% T20 = interp2(k0hd,hshd,T20Ref,k0.*h_d+fill,  h_s./h_d+fill,'makima',nan);



for ih = length(h_s):-1:1
    k0s(ih,:) = findWaveNumbers(w0,h_s(ih),0,0);
    for iw = length(w0):-1:1
        %% [1] linear waves coefficients and wavenumbers
        [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s(ih),w0(iw),0);
        
        %% [2] super-harmonic wave coefficients and wavenumbers
        [~,T_2m,~,~] = Free_waves_super_Harmonic(h_d,h_s(ih),k_nv(1),k_msv(1),R_n(1),T_m(1),0);
        
        R0(ih,iw) = abs(R_n(1));
        T0(ih,iw) = abs(T_m(1));
        T20(ih,iw) = abs(T_2m(1));
        
    end
end

sig_s = tanh(k0s.*h_s);
heta_b = k0s.*(3-sig_s.^2)./(4*sig_s.^3).*T0.^2;
heta_f = 2*w0.^2./g.*T20  ;
heta_ratio = heta_f./heta_b;
heta_ratio(heta_ratio>2)=nan;
gfz = figure('color','w');
contourf(h_s./h_d+0*f,1./f+0*h_s,heta_ratio,'ShowText','on')
xlabel('h_s/h_d [m]'), ylabel('T [s]'); colorbar
title('$$|\eta_{0Tf}^{(22,0)}|/|\eta_{0Tb}^{(22,0)}|$$, ($$h_d = 10$$\,m)','interpreter','latex','fontsize',14)


hfR = figure('color','w');
contourf(h_s./h_d+0*f,1./f+0*h_s,R0./T0,'ShowText','on')
xlabel('h_s/h_d [m]'), ylabel('T [s]'); colorbar
title('$$|\eta_{0R}^{(11,0)}|/|\eta_{0T}^{(11,0)}|$$, ($$h_d = 10$$\,m)','interpreter','latex','fontsize',14)


hfT = figure('color','w');hfT.Position(3) = 1600;
subplot(1,3,1),contourf(h_s./h_d+0*f,1./f+0*h_s,R0);colorbar('east','FontWeight','bold','Color','w')
set(gca,'XScale','log','YScale','log');xlabel('h_s/h_d [m]'), ylabel('T [s]');
title('$$|R_0|$$','interpreter','latex','fontsize',14)
subplot(1,3,2),contourf(h_s./h_d+0*f,1./f+0*h_s,T0);colorbar('east','FontWeight','bold','Color','w')
set(gca,'XScale','log','YScale','log');xlabel('h_s/h_d [m]'), ylabel('T [s]');
title('$$|T_0|$$','interpreter','latex','fontsize',14)
T20NaN = T20; T20NaN(T20NaN>25)=nan;
subplot(1,3,3),contourf(h_s./h_d+0*f,1./f+0*h_s,T20NaN),colorbar('east','FontWeight','bold','Color','w')
set(gca,'XScale','log','YScale','log');xlabel('h_s/h_d [m]'), ylabel('T [s]');
title('$$|T_{20}|$$','interpreter','latex','fontsize',14)


figure('color','w');
contourf(h_s./h_d+0*f,1./f+0*h_s,heta_f);colorbar;title('heta_f');
figure('color','w');
contourf(h_s./h_d+0*f,1./f+0*h_s,heta_b);colorbar;title('heta_b');


% hfT = figure('color','w');hfT.Position(3) = 1600;
% [X,Y] = meshgrid(linspace(.3,6,200),linspace(.2,1,200));
% subplot(1,3,1),contourf(X,Y,interp2(k0hd,hshd,R0Ref,X,Y,'makima'));colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|R_0|$$','interpreter','latex','fontsize',14)
% hold on; plot(k0hd,hshd,'ko')
% subplot(1,3,2),contourf(X,Y,interp2(k0hd,hshd,T0Ref,X,Y,'makima'));colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|T_0|$$','interpreter','latex','fontsize',14)
% hold on; plot(k0hd,hshd,'ko')
% subplot(1,3,3),contourf(X,Y,interp2(k0hd,hshd,T20Ref,X,Y,'makima')),colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|T_{20}|$$','interpreter','latex','fontsize',14)
% hold on; plot(k0hd,hshd,'ko')


% hfY = figure('color','w');hfY.Position(3) = 1600;
% [X,Y] = meshgrid(linspace(.3,6,200),linspace(.2,1,200));
% subplot(1,3,1),contourf(k0.*h_d+fill,  h_s./h_d+fill,abs(R0_Yan));colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|R_0|$$','interpreter','latex','fontsize',14)
% subplot(1,3,2),contourf(k0.*h_d+fill,  h_s./h_d+fill,abs(T0_Yan));colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|T_0|$$','interpreter','latex','fontsize',14)
% ha = subplot(1,3,3);contourf(k0.*h_d+fill,  h_s./h_d+fill,abs(T20_Yan),[20,10,5,3,1,0]);
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|T_{20}|$$','interpreter','latex','fontsize',14)
% hc = ha.Children(end);
% % hc.LevelList = [20,10,5,3,1,0];
% % colorbar('east','FontWeight','bold','Color','w')
% colorbar
% % hold on;plot(k0.*h_d+fill,  h_s./h_d+fill,'k.')

% contourf(k0hd,hshd , interp2(k0hd,hshd,T20Ref,k0.*h_d+fill,h_s./h_d+fill,'makima'))