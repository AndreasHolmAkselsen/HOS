clear
addpath .\step_functions_YanLi

h_s = logspace(log10(2),log10(9.5),50)';
f = logspace(log10(.1),log10(.36),50);
% h_s = (2:.5:9.5)';
% f = (.1:.01:.36);
h_d = 10;
g = 9.81;



NMode = 100; % hardly any difference from =0.

omega_0 = 2*pi*f;
k_0 = findWaveNumbers(omega_0,h_d,0,0) + 0*omega_0;
for i_hs = size(h_s,1):-1:1
    k_0s(i_hs,:) = findWaveNumbers(omega_0,h_s(i_hs),0,0);
end
k0hd = k_0.*h_d;

%% from Li's second_order_step_function.m
%-----------%-----------%
%      sub-harmonic     % incoming and reflected
%-----------%-----------%
c_g0              = 0.5.*omega_0./k_0.*(1+2.*k0hd./sinh(2.*k0hd));
term1B            = (2.*g.*h_d-c_g0.^2 )./2./sinh(2.*k0hd) + 2.*g.*c_g0./omega_0;
B_d               = -1 ./ (4.*(g.*h_d-c_g0.^2)) .*term1B;

%-----------%-----------%
%    sub--harmonic      % transmitted
%-----------%-----------%
k0hs              = k_0s.*h_s;
c_g0s             = 0.5.*omega_0./k_0s.*(1+2.*k0hs./sinh(2.*k0hs));
term1Bs            = (2.*g.*h_s-c_g0s.^2 )./2./sinh(2.*k0hs) + 2.*g.*c_g0s./omega_0;
B_s               = -1 ./ (4.*(g.*h_s-c_g0s.^2)) .*term1Bs;

%-----------%-----------%
%  sub-harmonic  free   %
%-----------%-----------%
for i_hs = size(h_s,1):-1:1
    for i_w = size(h_s,2):-1:1
%         [R_0(i_hs,i_w),T_0(i_hs,i_w)] = monochramonic_coefficient_final(h_d,h_s(i_hs,i_w),omega_0(1,i_w),0);
        [a,b] = monochramonic_coefficient_final(h_d,h_s(i_hs,i_w),omega_0(1,i_w),NMode);
        [R_0(i_hs,i_w),T_0(i_hs,i_w)] = deal(a(1),b(1));
    end
end
F_1               = k_0s.*B_s.*h_s.*g.*(abs(T_0)).^2./c_g0s - (1-(abs(R_0)).^2).*B_d.*h_d.*k_0*g./c_g0;  % from the momentum
F_2               = k_0s.*B_s.*(abs(T_0)).^2 - (1+(abs(R_0)).^2).*B_d.*k_0; 
B_Rf              = (F_1-sqrt(g*h_s).*F_2)./(-(sqrt(g*h_d)+sqrt(g*h_s)).*k_0);
B_Tf              = (F_1+sqrt(g*h_d).*F_2)./(-(sqrt(g*h_d)+sqrt(g*h_s)).*k_0s);



[H_S,F] = ndgrid(h_s ,f);

hfT = figure('color','w','name',['N modes: ',num2str(NMode)]);
contourf(H_S./h_d,1./F,abs(B_Tf)./abs(B_s),'ShowText','on')
xlabel('h_s/h_d'), ylabel('T [s]'); colorbar
title('$$|\eta_{T,f}^{(20,1)}|/|\eta_{T,b}^{(20,1)}|$$, ($$h_d = 10$$\,m)','interpreter','latex','fontsize',14)

hfR = figure('color','w','name',['N modes: ',num2str(NMode)]);
contourf(H_S./h_d,1./F,abs(B_Rf)./abs(B_d),'ShowText','on')
xlabel('h_s/h_d'), ylabel('T [s]'); colorbar
title('$$|\eta_{R,f}^{(20,1)}|/|\eta_{R,b}^{(20,1)}|$$, ($$h_d = 10$$\,m)','interpreter','latex','fontsize',14)

return
savefig(hfT,'.\memo\figures\zetaT201')
export_fig(hfT,'.\memo\figures\zetaT201','-png','-pdf','-m2')
savefig(hfR,'.\memo\figures\zetaR201')
export_fig(hfR,'.\memo\figures\zetaR201','-png','-pdf','-m2')
