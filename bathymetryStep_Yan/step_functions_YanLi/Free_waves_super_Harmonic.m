function [R_n,T_m,k_2nv,k_2msv] = Free_waves_super_Harmonic(h_d,h_s,k_0,k_0s,R_01st,T_01st,N_mode)
% --------------------------------- %
% by Yan Li, contact yan.li@ntnu.no %
% --------------------------------- %

% note here that now k_0 is 2 times the carrier wave number
%% finite water range: 0<< k h < pi
g      = 9.81;

%% return k^r_n which are the reflected waves in water depth h_d: solving omega^2 = - g k_n tan (k_n*h_d)
% return k_m which are the reflected waves in water depth h_s: solving omega^2 = - g k_m tan (k_m*h)
omega_0 = sqrt(g*k_0*tanh(k_0*h_d));
k_20     = 2*k_0;
omega    = 2*omega_0;

% const_1  = omega^2*h_d/g;
% const_2  = omega^2*h_s/g;
% cof_01   = 1;
% cof_02   = 1;
% if const_1/(0.75*pi) >0.9
%     cof_01   = 0.75;
%    if const_1/(0.56*pi) > tan(0.56*pi)
%        cof_01= 0.56;
%    end
% end
% 
% if const_2/(0.75*pi) >0.9
%        cof_02     = 0.75;
%    if  const_2/(0.56*pi) > tan(0.56*pi)
%        cof_02     =  0.56;
%    end
% end

% % AHA
k_d = findWaveNumbers(omega,h_d,0,N_mode-1);
k_2nv0 = k_d(1);
alpha_2nrvp = 1i*k_d(2:end).';
k_s = findWaveNumbers(omega,h_s,0,N_mode);
k_20s = k_s(1);
alpha_2mrvp = 1i*k_s(2:end).';

k_2nv        = [k_2nv0 1i*alpha_2nrvp];
k_20          = k_2nv0;
k_2mv         = 1i*alpha_2mrvp;
k_2msv        = [k_20s k_2mv];  % of size N_mode+1

if min([alpha_2nrvp k_2nv0 alpha_2mrvp k_20s])<0
    fprintf('Attention: find a negative wave number!!')
end
sgn = 1;
if k_20s*h_s>2*pi
    R_n(1) = 0;
    R_n(1:N_mode) = 0;
    T_m(1) = 1;
    T_m(2:N_mode+1) = 0;
    sgn =0;
end

if sgn 
k_DtH   = k_2nv*h_d;
k_StH   = k_2msv*h_s;
ch_D    = cosh(k_DtH);
ch_S    = cosh(k_StH);

%% 
[k_nmat, k_mmat] = meshgrid(k_2nv,k_2msv);
Int_FnFn         = (sinh(2*k_DtH)+2*k_DtH)./4./k_2nv./ch_D.^2 ; % size(k_nv)
Int_GmGm         = (sinh(2*k_StH)+2*k_StH)./4./k_2msv./ch_S.^2;
sh_1             = sinh(k_nmat*h_d+k_mmat*h_s)./(k_nmat+k_mmat);
sh_2             = sinh(k_nmat*h_d-k_mmat*h_s)./(k_nmat-k_mmat);
sh_3             = 2*k_nmat./(k_nmat.^2-k_mmat.^2).*sinh(k_nmat*(h_d-h_s));
Int_FnGm         = (sh_1+sh_2-sh_3)/2./cosh(k_nmat*h_d)./cosh(k_mmat*h_s);

Diag_Rn          = [(-1i*k_20*Int_FnFn(1)) Int_FnFn(2:end).*alpha_2nrvp];
% rhs_F0Fn(1)      = Diag_Rn(1);
% rhs_F0Fn(2:length(k_nv)) = 0; 

%% vector on the right hand side
%%%-----------------%%%
%-%  incoming waves %-%
%%%-----------------%%%
sh_4th_d        = (sinh(k_0*h_d)).^4;
coeff_sum_II    = 3*omega_0/8./sh_4th_d.*cosh(2.*k_0.*h_d);  % no amplitude is considered
Int_d_II_gm     = obtain_integral_FFGG(h_d,h_s,2*k_0,k_2msv,2);
Int_d_II_fn     = obtain_integral_FFGG(h_d,h_s,2*k_0,k_2nv,1);
rhs_d_gm        = coeff_sum_II.*Int_d_II_gm;
rhs_d_fn        = coeff_sum_II.*Int_d_II_fn;

%%%-----------------%%%
%-% reflected waves %-%
%%%-----------------%%%
sh_4th_d        = (sinh(k_0*h_d)).^4;
coeff_sum_RR    = 3*R_01st^2*omega_0/8./sh_4th_d.*cosh(2.*k_0.*h_d);  % no amplitude is considered
Int_d_RR_gm     = obtain_integral_FFGG(h_d,h_s,2*k_0,k_2msv,2);
Int_d_RR_fn     = obtain_integral_FFGG(h_d,h_s,2*k_0,k_2nv,1);
rhs_d_RR_gm     = coeff_sum_RR.*Int_d_RR_gm;
rhs_d_RR_fn     = coeff_sum_RR.*Int_d_RR_fn;

%%%-----------------%%%
%-%transmitted waves%-%
%%%-----------------%%%
sh_4th_s         = (sinh(k_0s*h_s)).^4;
coeff_sum_II     = 3 *T_01st^2* omega_0/8./sh_4th_s.*cosh(2.*k_0s.*h_s);    % no amplitude is considered
Int_s_II_gm      = obtain_integral_FFGG(h_d,h_s,2*k_0s,k_2msv,3);
Int_s_II_fn      = obtain_integral_FFGG(h_d,h_s,k_2nv,2*k_0s,2);
rhs_s_gm         = coeff_sum_II.*Int_s_II_gm;
rhs_s_fn         = coeff_sum_II.*Int_s_II_fn;
%% Solve Rn+Tm here: coefficents matrix of size (N+M)*(N+M) 
c_Rn             = [ Int_FnGm  (diag(-Int_GmGm)); diag(Diag_Rn) -1i*(k_mmat.').*(Int_FnGm.')]; 
Rhs_eq           = [(rhs_s_gm-rhs_d_gm-rhs_d_RR_gm).';...
                   (2*1i*k_0s*rhs_s_fn-2*1i*k_0*rhs_d_fn + 2*1i*k_0*rhs_d_RR_fn).'];
RT_mn            = lsqr(c_Rn,Rhs_eq,10^-8,500);
R_n              = RT_mn(1:length(k_2nv)).';
T_m              = RT_mn(length(k_2nv)+1:end).';

end

% a = 0;
end

