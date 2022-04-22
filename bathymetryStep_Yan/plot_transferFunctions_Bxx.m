clear
% k0hd = .25:.5:6;
% hs_hd = .5:.05:1';
[k0hd,hshd] = meshgrid( logspace(log10(.3),log10(6.0),40), logspace(log10(.2),log10(.99),45));

h_d = 10;
g = 9.81;

NMode = 100; % don't see any difference from 0.

addpath .\step_functions_YanLi

h_s = hshd.*h_d;
k_0 = k0hd./h_d;
omega_0 = sqrt(g.*k_0.*tanh(k_0.*h_d));

for i_hs = size(h_s,1):-1:1
    k_0s(i_hs,:) = findWaveNumbers(omega_0(1,:),h_s(i_hs,1),0,0);
end


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

%%%%%%%%%



funcs = {'B_Rf','B_d','B_Tf','B_s'};
names = {'$$B_R^f$$','$$B_d$$','$$B_T^f$$','$$B_s$$'};
contours = {
  [-2,-1,-.5,0]
  [-5,-3,-2,-1,-.5,0]
  [10,5,3,1,0]
  [-20,-10,-5,-3,-2,-1,-.5,0]
};
hfY = figure('color','w');hfY.Position(3:4) = [1200,1000];
for i = 1:length(funcs)
    subplot(2,2,i),contourf(k0hd,hshd,eval(funcs{i}),contours{i},'ShowText','on');colorbar('east','FontWeight','bold','Color','w')
    set(gca,'XScale','log','YScale','log');xlabel('$$k_0 h_d$$','interpreter','latex','fontsize',14); ylabel('$$h_s/h_d$$','interpreter','latex','fontsize',14); 
    title(names{i},'interpreter','latex','fontsize',14)
end











% B_Rf_2  = ( g.*k_0s.*h_s.*B_s.*abs(T_0).^2./c_g0s...
%          - g.*k_0.*h_d.*B_d.*(1-abs(R_0).^2)./c_g0...
%           -sqrt(g*h_s).*(k_0s.*B_s.*abs(T_0).^2  - k_0.*B_d.*(abs(R_0).^2+1)) ...
%                 )./(-(sqrt(g*h_d)+sqrt(g*h_s)).*k_0);
% 
% B_Rf_paper= ( g*k_0s.*h_s.*B_s.*abs(T_0).^2./c_g0s ...
%             - g*k_0*h_d.*B_d.*(1-abs(R_0).^2)./c_g0 ...
%             - sqrt(g*h_s).*( k_0s.*B_s.*abs(T_0).^2 - k_0.*B_d.*(abs(R_0).^2+1) )...
%            ) ./ ( ( -sqrt(g*h_d) - sqrt(g*h_s) ).*k_0);
%        
%  
% B_Tf_2    = ( -  g.*k_0s.*h_s.*B_s.*abs(T_0).^2./c_g0s...
%             + g.*k_0.*h_d.*B_d.*(1-abs(R_0).^2)./c_g0...
%             - sqrt(g*h_d).*(k_0s.*B_s.*abs(T_0).^2 - k_0.*B_d.*(abs(R_0).^2+1) )...
%           )./( (sqrt(g*h_d)+sqrt(g*h_s)).*k_0s);
%        
%         
% B_Tf_paper= ( g*k_0s.*h_s.*B_s.*abs(T_0).^2./c_g0s ...
%             + g*k_0.*h_d.*B_d.*(1-abs(R_0).^2)./c_g0 ...
%             + sqrt(g.*h_d).*( k_0s.*B_s.*abs(T_0).^2 - k_0.*B_d.*(abs(R_0).^2+1) )...
%            ) ./( ( sqrt(g.*h_d) + sqrt(g.*h_s) ).*k_0s);       
% 
% 
% B_Tf_code = ( k_0s.*h_s.*B_s.*abs(T_0).^2 ...
%             - k_0.*h_d.*B_d.*(1-abs(R_0).^2)...
%             + sqrt(g.*h_d).*( k_0s.*B_s.*abs(T_0).^2  - k_0.*B_d.*(abs(R_0).^2+1) )...
%            ) ./( ( -sqrt(g.*h_d) - sqrt(g.*h_s) ).*k_0s);
% 
% B_Rf_code = ( k_0s.*h_s.*B_s.*abs(T_0).^2 ...
%             - k_0.*h_d.*B_d.*(1-abs(R_0).^2)...
%             - sqrt(g.*h_s).*(k_0s.*B_s.*abs(T_0).^2-k_0.*B_d.*(1+abs(R_0).^2))...
%            )./( (-sqrt(g.*h_d)-sqrt(g.*h_s) ).*k_0);            
%             
% 
% % B_Tf_edit = -(  g*k_0s.*h_s.*B_s.*abs(T_0).^2./c_g0s ...
% %             - g*k_0*h_d.*B_d.*(1-abs(R_0).^2)./c_g0 ...
% %             + sqrt(g*h_d).*( k_0s.*B_s.*abs(T_0).^2 - k_0.*B_d.*(abs(R_0).^2+1) )...
% %            ) ./( ( sqrt(g*h_d) + sqrt(g*h_s) ).*k_0s);
%        