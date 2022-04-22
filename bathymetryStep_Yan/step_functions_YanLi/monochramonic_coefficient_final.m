function [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s,omega,N_mode)
% --------------------------------- %
% by Yan Li, contact yan.li@ntnu.no %
% --------------------------------- %
%%
% R_n: size of 1*N_mode for reflection
% T_m: Size of 1*N_mode+1 for transmission

%% finite water range: 0<< k h < pi
% omega  = 2*pi;
% h_d    = 0.3; 
% h_s    = 0.15;
g      = 9.81;
% N_mode = 32 ; % k_n < N_mode*pi

%% return k^r_n which are the reflected waves in water depth h_d: solving omega^2 = - g k_n tan (k_n*h_d)
% return k_m which are the reflected waves in water depth h_s: solving omega^2 = - g k_m tan (k_m*h)
% omega = sqrt(g*k_0*tanh(k_0*h_d));

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

% AHA
k_d = findWaveNumbers(omega,h_d,0,N_mode-1);
k_nv0 = k_d(1);
k_nrvp = 1i*k_d(2:end).';
k_s = findWaveNumbers(omega,h_s,0,N_mode-1);
k_0s = k_s(1);
k_mrvp = 1i*k_s(2:end).';
k_nv     = [k_nv0 1i*k_nrvp];
k_0      = k_nv0;
k_mv     = 1i*k_mrvp;
k_msv    = [k_0s k_mv];  % of size N_mode+1

if min([k_nrvp k_nv0 k_mrvp k_0s])<0
    fprintf('Attention: find a negative wave number!!')
end
sgn = 1;
if k_0s*h_s>2*pi
    R_n(1) = 0;
    R_n(1:N_mode) = 0;
    T_m(1) = 1;
    T_m(2:N_mode+1) = 0;
    sgn =0;
end

if sgn 
k_DtH   = k_nv*h_d;
k_StH   = k_msv*h_s;
ch_D    = cosh(k_DtH);
ch_S    = cosh(k_StH);

%% 
[k_nmat, k_mmat] = meshgrid(k_nv,k_msv);
Int_FnFn         = (sinh(2*k_DtH)+2*k_DtH)./4./k_nv./ch_D.^2 ; % size(k_nv)
Int_GmGm         = (sinh(2*k_StH)+2*k_StH)./4./k_msv./ch_S.^2;
sh_1             = sinh(k_nmat*h_d+k_mmat*h_s)./(k_nmat+k_mmat);
sh_2             = sinh(k_nmat*h_d-k_mmat*h_s)./(k_nmat-k_mmat);
sh_3             = 2*k_nmat./(k_nmat.^2-k_mmat.^2).*sinh(k_nmat*(h_d-h_s));
Int_FnGm         = (sh_1+sh_2-sh_3)/2./cosh(k_nmat*h_d)./cosh(k_mmat*h_s);

Diag_Rn          = [(-1i*k_0*Int_FnFn(1)) Int_FnFn(2:end).*k_nrvp];
rhs_F0Fn(1)      = Diag_Rn(1);
rhs_F0Fn(2:length(k_nv)) = 0; 

%% Solve Rn+Tm here: coefficents matrix of size (N+M)*(N+M) 
c_Rn             = [ Int_FnGm  (diag(-Int_GmGm)); diag(Diag_Rn) -1i*(k_mmat.').*(Int_FnGm.')]; 
Rhs_eq           = [(-Int_FnGm(:,1)); Diag_Rn(1);0.*k_nrvp'];
RT_mn            = lsqr(c_Rn,Rhs_eq,10^-8,500);
R_n              = RT_mn(1:length(k_nv)).';
T_m              = RT_mn(length(k_nv)+1:end).';
end

end

