function [int] = obtain_integral_FFGG(h_d,h_s,k_1,k_2,case_sign) 
% --------------------------------- %
% by Yan Li, contact yan.li@ntnu.no %
% --------------------------------- %
%%
% Note this function has input parameters that do not involve
% alpha_n/alpha_m, but should be alway k_n/k_m
% the sign is used to indicate the water depth appeared in the sinh/cosh functions only. 
switch case_sign
    % both deepwater
    case 1
        sh1 = sinh(h_d.*(k_1+k_2));
        sh2 = sinh(h_d.*(k_1-k_2));
        int = 1/2./cosh(k_1.*h_d)./cosh(k_2.*h_d)...
              .*(sh1./(k_1+k_2)+sh2./(k_1-k_2));
        
    % one deep+one shallow
    case 2
        
%         sh_1             = sinh(k_nmat*h_d+k_mmat*h_s)./(k_nmat+k_mmat);
%         sh_2             = sinh(k_nmat*h_d-k_mmat*h_s)./(k_nmat-k_mmat);
%         sh_3             = 2*k_nmat./(k_nmat.^2-k_mmat.^2).*sinh(k_nmat*(h_d-h_s));
%         Int_FnGm         = (sh_1+sh_2-sh_3)/2./cosh(k_nmat*h_d)./cosh(k_mmat*h_s);
        
        sh1 = sinh(k_1.*h_d+k_2.*h_s);
        sh2 = sinh(k_1.*h_d-k_2.*h_s);
        sh3 = sinh(k_1.*(h_d-h_s));
        int = 1/2./cosh(k_1.*h_d)./cosh(k_2.*h_s)...
              .*(sh1./(k_1+k_2)+sh2./(k_1-k_2)-2*k_1.*sh3./(k_1.^2-k_2.^2));
          
    % both shallow
    case 3 
        sh1 = sinh(h_s.*(k_1+k_2));
        sh2 = sinh(h_s.*(k_1-k_2));
        int = 1/2./cosh(k_1.*h_s)./cosh(k_2.*h_s)...
              .*(sh1./(k_1+k_2)+sh2./(k_1-k_2));
end
end

%% test below
% N  = 100000;
% zv = linspace(-h_s,0,N);
% dz = h_s/(N-1);
% [k_2mat,zv_mat] = meshgrid(k_2,zv);
% 
% int_test = dz*trapz(cosh(k_1.*((zv_mat+h_d))).*cosh(k_2mat.*((zv_mat+h_s))),1)...
%            ./cosh(k_1*h_d)./cosh(k_2.*h_s);

