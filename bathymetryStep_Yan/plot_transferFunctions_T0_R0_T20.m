clear
% k0hd = .25:.5:6;
% hs_hd = .5:.05:1';
[k0hd,hshd] = meshgrid( logspace(log10(.25),log10(6.0),40), logspace(log10(.2),log10(.99),45));

h_d = 10;
g = 9.81;

% NMode = 0;
NMode = 200;

addpath .\step_functions_YanLi

h_s = hshd.*h_d;
k0 = k0hd./h_d;
w0 = sqrt(g*k0.*tanh(k0*h_d));

for i = size(k0,1):-1:1
    for j = size(k0,2):-1:1
%         k0s(i,j) = findWaveNumbers(w0(i,j),h_s(i,j),0,0);
        
         %% [1] linear waves coefficients and wavenumbers
        [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s(i,j),w0(i,j),NMode);
        
        %% [2] super-harmonic wave coefficients and wavenumbers
        [~,T_2m,~,~] = Free_waves_super_Harmonic(h_d,h_s(i,j),k_nv(1),k_msv(1),R_n(1),T_m(1),NMode);
        
        R0(i,j) = R_n(1);
        T0(i,j) = T_m(1);
        T20(i,j) = T_2m(1);       
    end
end



hf = figure('color','w','name',['N modes: ',num2str(NMode)]);hf.Position(3) = 1600;
subplot(1,3,1),contourf(k0hd,hshd,abs(R0),[.3,.2,.15,.1,.05,0],'ShowText','on');colorbar('east','FontWeight','bold','Color','w')
set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
title('$$|R_0|$$','interpreter','latex','fontsize',14)
subplot(1,3,2),contourf(k0hd,hshd,abs(T0),[1.3,1.2,1.1,1.05,1,.95,0],'ShowText','on');colorbar('east','FontWeight','bold','Color','w')
set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
title('$$|T_0|$$','interpreter','latex','fontsize',14)
T20NaN = T20; T20NaN(abs(T20)>40)=nan;
ha = subplot(1,3,3);contourf(k0hd,hshd,abs(T20NaN),[20,10,5,3,1,0],'ShowText','on');
set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
title('$$|T_{20}|$$','interpreter','latex','fontsize',14)
hc = ha.Children(end);
colorbar('east','FontWeight','bold','Color','w')
% colorbar

return
savefig(hf,'.\memo\figures\T_contour_Yan')
export_fig(hf,'.\memo\figures\T_contour_Yan','-png','-pdf','-m2')
