%% Obtaining the evolution of a Gaussian packet over a step, correct to second order in wave steepness
% --------------------------------- %
% by Yan Li, contact yan.li@ntnu.no %
% --------------------------------- %
 
%% [0] return the inputs: for more explanations of the paramters, check the <myinputs> function
% [h_d,h_s,x_focus,t_focus,omega_0,delta_x,eps,N_f,t_step] = myinputs() ; 
cv          = [0.55 0.8 0.06 0 0.06] ;
h_d     = cv(1); 
h_s     = h_d -0.35; 
f_0     = cv(2);
omega_0 = 2*pi*f_0;

N_mode  = 32; 
fs_wg       = 128;
time        = 1/fs_wg:1/fs_wg:64;
[ipt]   = build_ipt_str(cv,time); % just to build up an input structure for the super-harmonic coefficients
% ipt.h_d     = cv(1);
% ipt.h_s     = cv(1) -0.35;
% ipt.f0      = cv(2);
% ipt.omega_0 = 2*pi*cv(2);
% ipt.delta   = cv(3);
% ipt.phase   = cv(4);
% ipt.epsilon = cv(5);
% ipt.t_vec   = t_vec;
% ipt.x_0     = 0.9;
% ipt.T_0     = 32 ;

% AHA
heaviside = @(x) (x>0);

%% [1] linear waves coefficients and wavenumbers 
[R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s,omega_0,32); 

%% [2] super-harmonic wave coefficients and wavenumbers
[R_2n,T_2m,k_2nv,k_2msv] = Free_waves_super_Harmonic(h_d,h_s,k_nv(1),k_msv(1),R_n(1),T_m(1),32);

%% [3] obtaining the elevation at different times and locations
[ipt]   = build_ipt_str(cv,time);
x_pvec = linspace(-5,5,256) ; % position vectors; the step is located at x =0 !!!

for ii = length(x_pvec):-1:1
    ii
    zeta = second_order_step_function(x_pvec(ii),R_2n,T_2m,k_2nv,k_2msv,R_n,T_m,k_msv,k_nv,ipt) ; % total surface elevation; up to second order in wave steepness
    zeta_d(:,ii)     = zeta.zeta_D;
    zeta_s(:,ii)     = zeta.zeta_S;
    zeta_1st(:,ii)   = zeta.zeta_1st;
    zeta_1stS(:,ii)  = zeta.zeta_1stS;
    %    zeta_Sp(:,i)  = zeta.zeta_S ;   % transmitted wavepacket
    zeta_superS(:,ii)= zeta.zeta_supS;
    zeta_superd(:,ii)= zeta.zeta_sup;
    
    zeta_subS(:,ii)  = zeta.zeta_subS;
    zeta_subd(:,ii)  = zeta.zeta_sub;
    zeta_errD(:,ii)  = zeta.zeta_errD;
end
delta_d    = heaviside(-x_pvec);
delta_s    = heaviside(x_pvec); 
zeta_all      = zeta_d.*delta_d      + zeta_s.*delta_s; % total elevation
zeta_allsuper = zeta_superd.*delta_d + zeta_superS.*delta_s; % super-harmonic
zeta_allsub   = zeta_subd.*delta_d   + zeta_subS.*delta_s; % sub-harmonic
zeta_1st_all0  = zeta_1st.*delta_d    + zeta_1stS.*delta_s; % linear 

%% [] plot

nx = 64 ;

close all
subplot(4,1,1)
plot(time*f_0,zeta_all(:,nx))
ylabel('$\zeta^{(1)}+\zeta^{(20)}+\zeta^{(22)} $', 'Interpreter','latex')

subplot(4,1,2)
plot(time*f_0,zeta_1st_all0(:,nx))
ylabel('$\zeta^{(1)} $', 'Interpreter','latex')

subplot(4,1,3)
plot(time*f_0,zeta_allsuper(:,nx))
ylabel('$\zeta^{(22)} $', 'Interpreter','latex')

subplot(4,1,4)
plot(time*f_0,zeta_allsub(:,nx))
xlabel('$t/T_p []$', 'Interpreter','latex')
ylabel('$\zeta^{(20)} $', 'Interpreter','latex')
