clear
addpath .\step_functions_YanLi

% h_s = (2:.5:9.5)';
% f = (.1:.01:.36);
% h_d = 10;
g = 9.81;

h_d = 10;%6;


h_s = logspace(log10(.2*h_d),log10(.95*h_d),45)';
f = logspace(log10(.1),log10(.36),45);


% NMode = 0;
NMode = 100;

w0 = 2*pi*f;
k0 = findWaveNumbers(w0,h_d,0,0) + 0*w0;
% k0s = 0*w0.*h_s;
% for i = 1:length(h_s)
%     k0s(i,:) = findWaveNumbers(w0,h_s(i),0,0);
% end

for ih = length(h_s):-1:1
    k0s(ih,:) = findWaveNumbers(w0,h_s(ih),0,0);
    for iw = length(f):-1:1
        %% [1] linear waves coefficients and wavenumbers
        [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s(ih),w0(iw),NMode);
        
        %% [2] super-harmonic wave coefficients and wavenumbers
        [~,T_2m,~,~] = Free_waves_super_Harmonic(h_d,h_s(ih),k_nv(1),k_msv(1),R_n(1),T_m(1),NMode);
        
        R0(ih,iw) = abs(R_n(1));
        T0(ih,iw) = abs(T_m(1));
        T20(ih,iw) = abs(T_2m(1));
        
    end
end

[H_S,F] = ndgrid(h_s ,f);

sig_s = tanh(k0s.*h_s);
heta_b = k0s.*(3-sig_s.^2)./(4*sig_s.^3).*T0.^2;
heta_f = 2*w0.^2./g.*T20  ;
heta_ratio = heta_f./heta_b;
heta_ratio(heta_ratio>2) = nan;

gfz = figure('color','w','name',['N modes: ',num2str(NMode)]);
contourf(H_S./h_d,1./F,heta_ratio)
xlabel('h_s/h_d'), ylabel('T [s]'); 
colorbar;%('east','FontWeight','bold','Color','w')
title(['$$|\zeta_{0Tf}^{(22,0)}|/|\zeta_{0Tb}^{(22,0)}|$$, ($$h_d = ',num2str(h_d),'$$\,m)'],'interpreter','latex','fontsize',14)


hfR = figure('color','w','name',['N modes: ',num2str(NMode)]);
contourf(H_S./h_d,1./F,R0./T0)
xlabel('h_s/h_d'), ylabel('T [s]'); 
colorbar;%('east','FontWeight','bold','Color','w')
title(['$$|\zeta_{0R}^{(11,0)}|/|\zeta_{0T}^{(11,0)}|$$, ($$h_d = ',num2str(h_d),'$$\,m)'],'interpreter','latex','fontsize',14)

return
post = ['_hd',num2str(h_d)]; post(post=='.')='p';
savefig(gfz,['.\memo\figures\zeta_f__b_Yan',post])
export_fig(gfz,['.\memo\figures\zeta_f__b_Yan',post],'-png','-pdf','-m2')
savefig(hfR,['.\memo\figures\zeta_R_Yan',post])
export_fig(hfR,['.\memo\figures\zeta_R_Yan',post],'-png','-pdf','-m2')

% hfT = figure('color','w');hfT.Position(3) = 1600;
% [X,Y] = meshgrid(linspace(.3,6,200),linspace(.2,1,200));
% fill = zeros(length(h_s),1);
% subplot(1,3,1),contourf(k0.*h_d+fill,  H_S./h_d,abs(R0));colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|R_0|$$','interpreter','latex','fontsize',14)
% subplot(1,3,2),contourf(k0.*h_d+fill,  H_S./h_d,abs(T0));colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|T_0|$$','interpreter','latex','fontsize',14)
% T20NaN = T20; T20NaN(abs(T20)>25)=nan;
% ha = subplot(1,3,3);contourf(k0.*h_d+fill,  H_S./h_d,abs(T20NaN));%,[20,10,5,3,1,0]);
% set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
% title('$$|T_{20}|$$','interpreter','latex','fontsize',14)
% hc = ha.Children(end);
% % colorbar('east','FontWeight','bold','Color','w')
% colorbar
% 
% 
% hfT = figure('color','w');hfT.Position(3) = 1600;
% subplot(1,3,1),contourf(H_S./h_d,1./F,R0);colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('h_s/h_d [m]'), ylabel('T [s]');
% title('$$|R_0|$$','interpreter','latex','fontsize',14)
% subplot(1,3,2),contourf(H_S./h_d,1./F,T0);colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('h_s/h_d [m]'), ylabel('T [s]');
% title('$$|T_0|$$','interpreter','latex','fontsize',14)
% T20NaN = T20; T20NaN(T20NaN>5)=nan;
% subplot(1,3,3),contourf(H_S./h_d,1./F,T20NaN),colorbar('east','FontWeight','bold','Color','w')
% set(gca,'XScale','log','YScale','log');xlabel('h_s/h_d [m]'), ylabel('T [s]');
% title('$$|T_{20}|$$','interpreter','latex','fontsize',14)
% 
% gfz = figure('color','w');
% contourf(H_S./h_d+0*F,1./F,heta_f);colorbar;title('heta_f');
% gfz = figure('color','w');
% contourf(H_S./h_d+0*F,1./F,heta_b);colorbar;title('heta_b');


