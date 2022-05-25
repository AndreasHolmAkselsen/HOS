function [zeta] = second_order_step_function(x_p,R_2n,T_2m,k_2nv,k_2msv,...
                                      R_n,T_m,k_msv,k_nv,phi_shift)
g        = 9.81;

% or use a function to place the inputs
[h_d,h_s,x_0,T_0,omega_0,delta_x,eps,N_f,t_step]  = myinputs();

k_0       = k_nv(1);
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
% S(omega)
% H_s         = 4*10^-3;
% d_omega     = 2*pi*f_s/N;
% % S_omega     = (H_s/4)^2/sigma/sqrt(2*pi)*exp(-(omega_vec -0).^2/2/sigma^2);
% % A0_p        = sqrt((H_s/4)^2/sigma/sqrt(2*pi)*d_omega*2);
% S_p         = (H_s/4)^2/sigma_omg/sqrt(2*pi);
% S_omega     = S_p*exp(-(omega_vec-omega_0).^2/2/sigma_omg^2);
% A0_p        = sqrt(2*sum(S_omega)*d_omega);  % coff = S_p/max(psd_zeta)
% coeff       = sqrt(24.3698);
% A0_p        = coeff*2*sqrt(2)*sqrt(pi)*A0_p*sigma_omg;
%-----------%-----------%
%     First order       %
%-----------%-----------%
% T_vec       =  (0:(N-1))*T_s;                % the slow time
% t_vec       = T_vec*c_g0/sqrt(2)/sigma_omg;  % the fast time
% 
% delta       = sqrt(2)*sigma_omg/c_g0;
% dt          = T_s*c_g0/sqrt(2)/sigma_omg;
% t_sum       = N*dt;

%% the fast and slow scales
t_vec       = (0:(N_f-1))*t_step;     % the slow time
T_vec       = t_vec/sigma_x;  % the fast time
X_p         = (x_p)/sigma_x;   % the slow X,the center is located at (-x_0)
X_0         = (x_0)/sigma_x;
x_0         = 0;
% T_0         = t_0;
T_0         = delta_x*T_0;
t_0         = 0;

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
B_d               = -1 ./ (4*k0hd.*(1-tanh(k0hd)/4/k0hd.*(1+2*k0hd./sinh(2*k0hd))^2)); 
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
B_s               = -1 ./ (4*k_0*h_s.*(1-tanh(k0hs)/4/k0hs.*(1+2*k0hs./sinh(2*k0hs))^2)); 
zeta_sub_T        =  B_s.*k_0.*A_envelopS.^2;

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

% plot(t_vec,zeta_sub+zeta_sup,'-k',t_vec,zeta_subS+zeta_supS,'--r',...
%      t_vec,zeta_2T_f,'-.b',t_vec,zeta_2R_f,'-.m') 
% set(gca,'xlim',[0 50])

%% now the second-order sub-harmonic FREE waves
%-----------%-----------%
%  sub-harmonic  free   %
%-----------%-----------%
c_gsubR           = sqrt(g*h_d);
c_gsubT           = sqrt(g*h_s);
gamma_h           = h_s/h_d;
F_1               = -(1+(abs(R_n(1)))^2)*B_d+(abs(T_m(1)))^2*B_s;
F_2               = -(1-(abs(R_n(1)))^2)*B_d+(abs(T_m(1)))^2*B_s*h_s/h_d;
B_Rf              = (sqrt(gamma_h)*F_1-c_gsubR/c_g0*F_2)/(1+sqrt(gamma_h));
B_Tf              = -(F_1+c_gsubR/c_g0*F_2)/(1+sqrt(gamma_h));
A_subRf           = A0_p*exp(-(-c_g0/c_gsubR*X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);
A_subTf           = A0_p*exp(-(c_g0/c_gsubT*X_p-X_0-c_g0.*(T_vec-T_0)).^2/2);
zeta_sub_Rf       = B_Rf*k_0.*A_subRf.^2;
zeta_sub_Tf       = B_Tf*k_0.*A_subTf.^2;

zeta.zeta_sub_Tf = zeta_sub_Tf; 
zeta.zeta_sub_Rf = zeta_sub_Rf;
%% evanescent waves
% zeta_Es_1         = sum(T_2m(2:end).*exp(1i*(k_2msv(2:end)*x_p - omega_0*t_vec)...
%                         + 1i*phi_shift),2);
% zeta_Ed_1         = sum(R_2n(2:end).*exp(1i*(-k_2nv(2:end)*x_p) + phi_shift),2);
% 
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

%% Surface elevation
zeta.zeta_1st   = zeta_1st + zeta_1stR;
zeta.zeta_sub   = zeta_sub_R + zeta_sub_I + zeta_sub_Rf;
zeta.zeta_sup   = zeta_sup_b + zeta_supR_b + zeta_2R_f;
zeta.zeta_D     = zeta.zeta_1st + zeta.zeta_sub+zeta.zeta_sup;
%
% zeta.zeta_S   = zeta_S;
zeta.zeta_1stS  = zeta_1stS ;
zeta.zeta_subS  = zeta_sub_T + zeta_sub_Tf;
zeta.zeta_supS  = zeta_supS  + zeta_2T_f ;
zeta.zeta_S     = zeta.zeta_1stS+ zeta.zeta_subS + zeta.zeta_supS;
zeta.zeta_errD  = zeta_errD;


end


