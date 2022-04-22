clear
% addpath .\step_functions_YanLi

nEv = 50;

k0shs = linspace(0.005,1.2,50);
% k0hs = .005:.01:1.2;
h_d = .8;
h_s = .15;

g = 9.81;

% k0 = k0hs./h_s;
% w0 = sqrt(g*k0.*tanh(k0*h_d));

k_0s = k0shs/h_s;
w0 = sqrt(g*k_0s.*tanh(k_0s*h_s));



for j = size(k_0s,2):-1:1
    %         k0s(i,j) = findWaveNumbers(w0(i,j),h_s(i,j),0,0);
    
    %% [1] linear waves coefficients and wavenumbers
    [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s,w0(1,j),nEv+1);
    
    [k_0(j),k_0s_check(j)] = deal(k_nv(1),k_msv(1));
    R0(j) = R_n(1);
    T0(j) = T_m(1);
end
fprintf('difference in k0s: %g\n',max(abs(k_0s-k_0s_check)))

 
hfk = figure('color','w');
% hpk = plot(k0shs.',R0.',k0shs.',T0.','LineWidth',1);
% xlabel('h_s k_{0s}'); ylabel('R_0');
hpk = plot(w0*sqrt(h_s/g),abs(R0),w0*sqrt(h_s/g),abs(T0),'LineWidth',1);
xlabel('\omega (h_s/g)^{1/2}'); ylabel('R_0');
grid on;
legend("h_s = "+h_s+"m",'location','northwest')
title("h_d = "+h_d+"m"); 
xlim([0,1.2])
ylim([0,1.4])
savefig('./benchmark_Gurerrez_w')
export_fig('./benchmark_Gurerrez_w','-png','-pdf')
