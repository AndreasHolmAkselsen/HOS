function [zeta] = second_order_step_function(x_p,R_2n,T_2m,k_2nv,k_2msv,...
                                      R_n,T_m,k_msv,k_nv,ipt)
% --------------------------------- %
% by Yan Li, contact yan.li@ntnu.no %
% --------------------------------- %
g        = 9.81;

% or use a function to place the inputs
h_d    = ipt.h_d;
h_s    = ipt.h_s;
x_0    = ipt.x_0;
T_0    = ipt.T_0;
omega_0= ipt.omega_0;
delta  = ipt.delta;
eps    = ipt.epsilon;
t_vec  = ipt.t_vec;
phi_shift = ipt.phase;


%%
k_0       = k_nv(1);
sigma     = 1/delta/k_0;
delta_x   = 1/sigma; 
sigma_x   = 1/delta_x;

%%
c_g0      = 0.5*omega_0/k_0*(1+2*k_0*h_d/sinh(2*k_0*h_d));

%-----------% --------- %
%   solve free waves    % 
%-----------% --------- %
k_0s        = k_msv(1);
c_g0s       = 0.5*omega_0/k_0s*(1+2*k_0s*h_s/sinh(2*k_0s*h_s));
% T_0         = 0;
%%  FFT  here 
% t_p       = 2*pi/omega_0;
% f_s       = 128;                  % sampling frequency (hz)
% t_step    = t_p/f_s;               % sampling period (s)
        
%% the PSD input

%% the fast and slow scales
% t_vec       = (0:(N_f-1))*t_step;     % the slow time
T_vec       = t_vec/sigma_x;  % the fast time
X_p         = (x_p)/sigma_x;   % the slow X,the center is located at (-x_0)
X_0         = (x_0)/sigma_x;
x_0         = 0;
t_0         = T_0;
T_0         = delta_x*t_0;


% AHA
heaviside = @(x) (x>0);

%%  
%-----------%-----------%
%  First-order incoming %
%-----------%-----------%
A0_p        = eps/k_0;
A_envelop   = A0_p*exp(-(X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);
zeta_1st    = A_envelop.*cos((k_0*(x_p-x_0)-omega_0*(t_vec-t_0))+phi_shift);

%-----------%-----------%
% First-order reflected %
%-----------%-----------%
phi_R          = atan(imag(R_n(1))/real(R_n(1))) +...
                  heaviside(-real(R_n(1))).*sign(imag(R_n(1)))*pi;
A_envelopR     = abs(R_n(1))*A0_p*exp(-(-X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);
zeta_1stR      = A_envelopR.*cos((-k_0*(x_p-x_0) - omega_0*(t_vec-t_0) + phi_R + phi_shift));
zeta.zeta_1stR = zeta_1stR;

%%
%-----------%-----------%
%First-order evanescentD%
%-----------%-----------%
% k_eD          = k_nv(2:end);
% c_g0_ED       = 0.5*omega_0./k_eD.*(1+2*k_eD*h_d./sinh(2*k_eD*h_d));
% c_g0_ED_mat   = repmat(c_g0_ED.',1,length(T_vec));
% T_vec_mat     = repmat(T_vec,length(c_g0_ED),1);
% A_envelop_ED  = A0_p*exp(-(c_g0./c_g0_ED_mat*X_p-X_0-c_g0.*(T_vec_mat-T_0)).^2/2);
% zeta_1st_ED   = sum(real(A_envelop_ED.*(exp(-1i*k_eD*x_p)).'.*cos(-omega_0*t_vec)));

%%
%-----------%-----------%
%    superharmonic      % BOUND: incoming and reflected
%-----------%-----------%
cos_sup     = cos(2*(k_0*(x_p-x_0)-omega_0*(t_vec-t_0))+2*phi_shift);
zeta_sup_b  = k_0*A_envelop.^2*cosh(k_0*h_d)*(2*(cosh(k_0*h_d))^2+1)/...
               4/(sinh(k_0*h_d))^3.*cos_sup;

phi_R       = atan(imag(R_n(1))/real(R_n(1)))+...
                heaviside(-real(R_n(1))).*sign(imag(R_n(1)))*pi;
A_envelopR  = abs(R_n(1))*A0_p*exp(-(-X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);
zeta_supR_b = k_0*A_envelopR.^2*cosh(k_0*h_d)*(2*(cosh(k_0*h_d))^2+1)/...
                4/(sinh(k_0*h_d))^3.*cos(2*(-k_0*(x_p-x_0) - omega_0*(t_vec-t_0) + phi_R)+2*phi_shift);

%-----------%-----------%
%      sub-harmonic     % incoming and reflected
%-----------%-----------%
%%
k0hd              = k_0*h_d;
c_g0              = 0.5*omega_0./k_0.*(1+2*k0hd./sinh(2*k0hd));
term1B            = (2*g*h_d-c_g0.^2 )/2./sinh(2*k0hd) + 2*g.*c_g0./omega_0;
B_d               = -1 ./ (4*(g*h_d-c_g0.^2)) .*term1B; 
zeta_sub_I        = B_d.*k_0.*A_envelop.^2;
zeta_sub_R        = B_d.*k_0.*A_envelopR.^2;

c_gd_err          = sqrt(g*h_d);           
A_err             = A0_p*exp(-(X_p-X_0-c_gd_err.*(T_vec-T_0)).^2/2);
zeta_errD         = -B_d.*k_0.*A_err.^2;

%% proceed to shallower linear and second-order bound waves
%-----------%-----------%
%        First          %
%-----------%-----------%
phi_T          = atan(imag(T_m(1))/real(T_m(1)))+...
                heaviside(-real(T_m(1))).*sign(imag(T_m(1)))*pi;
A_envelopS     = A0_p*abs(T_m(1))*exp(-(c_g0/c_g0s*X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);

zeta_1stS      = A_envelopS.*cos((k_0s*(x_p-x_0) - omega_0*(t_vec-t_0) + phi_T) + phi_shift);

%-----------%-----------%
%    superharmonic      %
%-----------%-----------%
cos_supS        = cos(2*(k_0s*(x_p-x_0)-omega_0*(t_vec-t_0) + phi_T) +  2*phi_shift);
zeta_supS       = k_0s*A_envelopS.^2*cosh(k_0s*h_s)*(2*(cosh(k_0s*h_s))^2+1)/4/...
                  (sinh(k_0s*h_s))^3.*cos_supS;

%-----------%-----------%
%    sub--harmonic      % transmitted
%-----------%-----------%
k0hs              = k_0s*h_s;

c_g0s             = 0.5*omega_0./k_0s.*(1+2*k0hs./sinh(2*k0hs));
term1B            = (2*g*h_s-c_g0s.^2 )/2./sinh(2*k0hs) + 2*g.*c_g0s./omega_0;
B_s               = -1 ./ (4*(g*h_s-c_g0s.^2)) .*term1B;
zeta_sub_T        =  B_s.*k_0s.*A_envelopS.^2;

%% evanescent waves 1st
k_eD              = k_nv(2:end);
c_g0_ED           = 0.5*omega_0./k_eD.*(1+2*k_eD*h_d./sinh(2*k_eD*h_d));
[T_matd,cgd_mat]  = meshgrid(T_vec,c_g0_ED);
k_eS              = k_msv(2:end);
c_g0_Es           = 0.5*omega_0./k_eS.*(1+2*k_eS*h_s./sinh(2*k_eS*h_s));
[T_mats,cgs_mat]   = meshgrid(T_vec,c_g0_Es);
A_ED              = A0_p*exp(-(-0.*c_g0./cgd_mat*X_p-X_0-c_g0.*(T_matd-T_0)).^2/2);
A_ES              = A0_p*exp(-(0.*c_g0./cgs_mat*X_p-X_0-c_g0.*(T_mats-T_0)).^2/2);
zeta_Ed_1         = sum(A_ED.*R_n(2:end).'.*exp(-1i*(k_nv(2:end).'*(x_p-x_0))...
                        + 1i*phi_shift));
zeta_Ed_1         = real(zeta_Ed_1.*exp(-1i*omega_0*(t_vec-t_0)));
zeta_Es_1         = real(sum(A_ES.*T_m(2:end).'.*exp(1i*(k_msv(2:end).'*(x_p-x_0))...
                        + 1i*phi_shift)).*exp(-1i*omega_0*(t_vec-t_0)));
                    
%% evanescent waves 2nd
if x_p == -0.1
    a = 0;
end
k_2eD              = k_2nv(2:end);
c_g0_2ED           = 0.5*2*omega_0./k_2eD.*(1+2*k_2eD*h_d./sinh(2*k_2eD*h_d));
[T_matd,cg2d_mat]  = meshgrid(T_vec,c_g0_2ED);
k_2eS              = k_2msv(2:end);
c_g0_2Es           = 0.5*2*omega_0./k_2eS.*(1+2*k_2eS*h_s./sinh(2*k_2eS*h_s));
[T_mats,cg2s_mat]   = meshgrid(T_vec,c_g0_2Es);
A_2ED              = A0_p^2*exp(-(-c_g0./cg2d_mat*X_p-X_0-c_g0.*(T_matd-T_0)).^2);
A_2ES              = A0_p^2*exp(-(c_g0./cg2s_mat*X_p-X_0-c_g0.*(T_mats-T_0)).^2);
zeta_2Ed_1         = sum(A_2ED.*R_2n(2:end).'.*exp(-1i*(k_2nv(2:end).'*(x_p-x_0))...
                        + 2*1i*phi_shift));
zeta_2Ed_1         = real(2*omega_0/g*zeta_2Ed_1.*exp(-1i*2*omega_0*(t_vec-t_0)));
zeta_2Es_1         = real(2*omega_0/g*sum(A_2ES.*T_2m(2:end).'.*exp(1i*(k_2msv(2:end).'*(x_p-x_0))...
                        + 2*1i*phi_shift)).*exp(-1i*2*omega_0*(t_vec-t_0)));

% zeta_2Es          = sum(real(...
%                            (2*1i*omega_0/9.81)*AT_sq_s(i,:).*zeta_Es_1.*...
%                            exp(-2*1i*omega_0.*t_vec).*...
%                            exp(-1i*2*phi_NB_s(i,:))...
%                            ));   
% 
% zeta_2Ed     = sum(real((2*1i*omega_0/9.81)*AT_sq_s.*zeta_Ed_1.*...
%                             exp(-2*1i*omega_0.*t_vec).*...
%                             exp(-2*1i*phi_NB_d)...
%                            ));

%% proceed to the second-order super-harmonic FREE waves
%-----------%-----------%
% superharmonic  free   %
%-----------%-----------%
phi_2T_f       =  atan(imag(T_2m(1))/real(T_2m(1)))+...
                  heaviside(-real(T_2m(1))).*sign(imag(T_2m(1)))*pi;
k_20s          =  k_2msv(1);
c_g20s         =  omega_0/k_20s*(1+2*k_20s*h_s/sinh(2*k_20s*h_s));
A_slowT        =  A0_p*exp(-(c_g0/c_g20s*(X_p)-X_0-c_g0.*(T_vec-T_0)).^2/2);
zeta_2T_f      =  2*omega_0/g*abs(T_2m(1)).*A_slowT.^2.*...
                  cos(k_20s.*(x_p-x_0)-2*omega_0*(t_vec-t_0)+phi_2T_f + 2*phi_shift);
zeta.zeta_2T_f =  zeta_2T_f;

phi_2R_f       =  atan(imag(R_2n(1))/real(R_2n(1)))+...
                  heaviside(-real(R_2n(1))).*sign(imag(R_2n(1)))*pi;
phi_2R_f(isnan(phi_2R_f))=0;

k_20d          =  k_2nv(1);
c_g20d         =  omega_0/k_20d*(1+2*k_20d*h_d/sinh(2*k_20d*h_d));

A_slowR        =  A0_p*exp(-(-c_g0/c_g20d*X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);
zeta_2R_f      =  2*omega_0/g*abs(R_2n(1)).*A_slowR.^2.*...
                  cos(-k_20d.*(x_p-x_0)-2*omega_0*(t_vec-t_0)+phi_2R_f  + 2*phi_shift);
zeta.zeta_2R_f =  zeta_2R_f;
%
zeta.T_vec     =  T_vec;
zeta.t_vec     =  t_vec;

%% now the second-order sub-harmonic FREE waves
%-----------%-----------%
%  sub-harmonic  free   %
%-----------%-----------%
c_gsubR           = sqrt(g*h_d);
c_gsubT           = sqrt(g*h_s);
F_1               = -(1-(abs(R_n(1)))^2)*k_0*h_d*B_d+(abs(T_m(1)))^2*B_s*k_0s*h_s;
F_2               = -(1+(abs(R_n(1)))^2)*B_d*k_0+(abs(T_m(1)))^2*B_s*k_0s;
B_Rf              = (F_1-sqrt(g*h_s).*F_2)./(-(sqrt(g*h_d)+sqrt(g*h_s)).*k_0);
B_Tf              = (F_1+sqrt(g*h_d).*F_2)./(-(sqrt(g*h_d)+sqrt(g*h_s)).*k_0s);
A_subRf           = A0_p*exp(-(-c_g0/c_gsubR*X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);
A_subTf           = A0_p*exp(-(c_g0/c_gsubT*X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);


zeta_sub_Rf       = B_Rf*k_0.*A_subRf.^2;
zeta_sub_Tf       = B_Tf*k_0s.*A_subTf.^2;

zeta.zeta_sub_Tf = zeta_sub_Tf; 
zeta.zeta_sub_Rf = zeta_sub_Rf;
%% evanescent waves
% no need to include as they are small

%% Surface elevation

zeta_Ed_1(isnan(zeta_Ed_1)) = 0;
zeta_Ed_1(isinf(zeta_Ed_1)) = 0;
zeta_Es_1(isnan(zeta_Es_1)) = 0;
zeta_Es_1(isinf(zeta_Es_1)) = 0;
zeta_2Ed_1(isnan(zeta_2Ed_1)) = 0;
zeta_2Ed_1(isinf(zeta_2Ed_1)) = 0;
zeta_2Es_1(isnan(zeta_2Es_1)) = 0;
zeta_2Es_1(isinf(zeta_2Es_1)) = 0;
if x_p == 0.1
    a = 0;
end
zeta.zeta_1st   = zeta_1st + zeta_1stR + zeta_Ed_1;
zeta.zeta_sub   = zeta_sub_R + zeta_sub_I + zeta_sub_Rf;
zeta.zeta_sup   = zeta_sup_b + zeta_supR_b + zeta_2R_f+ zeta_2Ed_1;
zeta.zeta_D     = zeta.zeta_1st + zeta.zeta_sub+zeta.zeta_sup;
zeta.zeta_1stS  = zeta_1stS + zeta_Es_1 ;
zeta.zeta_subS  = zeta_sub_T + zeta_sub_Tf;
zeta.zeta_supS  = zeta_supS  + zeta_2T_f + zeta_2Es_1;
zeta.zeta_S     = zeta.zeta_1stS+ zeta.zeta_subS + zeta.zeta_supS;
zeta.zeta_errD  = zeta_errD;
zeta.B_Tf       = B_Tf;
zeta.B_s        = B_s;
end


