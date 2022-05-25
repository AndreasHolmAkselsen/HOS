addpath O2StepTheory

h_d = 1;
h_s = 0.35;
T = 2.2987;
% T = 4.1733;
NMode = 100;


w = 2*pi/T;

k_d = findWaveNumbers(w,h_d,0,0)
k_s = findWaveNumbers(w,h_s,0,0)
lam_d = 2*pi/k_d
g = 9.81;

[R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s,w,NMode);
[~,T_2m,~,~] = Free_waves_super_Harmonic(h_d,h_s,k_nv(1),k_msv(1),R_n(1),T_m(1),NMode);



A = 0.05;
aRO1 = abs(R_n(1))*A
aO2 = A^2* 2*w/g*abs(T_2m(1))

abs(T_2m(1))



a2HOS = .020
T2HOS = a2HOS/(A^2*2*w/g)

% 
% A0_p        = eps/k_0;
% A_slowT        =  A0_p*exp(-(c_g0/c_g20s*(X_p)-X_0-c_g0.*(T_vec-T_0)).^2/2);
% zeta_2T_f      =  2*omega_0/g*abs(T_2m(1)).*A_slowT.^2.*...
%                   cos(k_20s.*(x_p-x_0)-2*omega_0*(t_vec-t_0)+phi_2T_f + 2*phi_shift);
% 
% zeta_sup_b  = k_0*A_envelop.^2*cosh(k_0*h_d)*(2*(cosh(k_0*h_d))^2+1)/...
%                4/(sinh(k_0*h_d))^3.*cos_sup;              