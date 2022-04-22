function [h_d,h_s,x_focus,t_focus,omega_0,delta_x,eps,N_f,t_step] = myinputs()
%% Time associated paramters
T      = 64;     % [s] time length
f_s    = 128;    % sampling frequency [Hz]
t_step = 1/f_s ; % time interval
N_f    = T/t_step+1; % total length of the time vector

%%

h_d        = 0.55;  % water depth on the deeper side
eps        = 0.043; % eps = k_0 A_p  with k_0 and A_p denoting hte characteristic wavelength and amplitude
delta      = 0.1;   % non-dimensional bandwidth parameter: delta=1/(k_0*sigma) where sigma is the standard deviation of a Gaussian evenlope
                    % length scale of the envelope is about 4*lambda_0 /delta_x
f_0        = 0.7;   % carrier wave frequency; [Hz]
x_focus    = 0.5;  %the position at focus
t_focus    = 32;   % the time at focus

%%
g          = 9.81;     
h_s        = h_d - 0.35; % water depth on the shallower deeper side
omega_0    = 2*pi*f_0; %[rad/s]
const_1    = omega_0^2*h_d/g;
k_0        = fsolve(@(x) -const_1./x+tanh(x),const_1)/h_d; % sovling the dispersion with a given f_0 and h_d
c_g0      = omega_0/2/k_0*(1+2*k_0*h_d/sinh(2*k_0*h_d));  % group speed
% lambda_0  = 2*pi/k_0;     % the carrier wavelength
% t_p       = 1/f_0;        % carrier wave period

% [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s,omega_0,32);

%  nondimensional parameters

%% Generation the inputs of a Gaussian wavepacket
% sigma=1/(k_0*delta);
% cg_0=g*(k_0*d*sech(k_0*d)^2+tanh(k_0*d))/(2*omega_0);
% a_0=epsilon/k_0;%amplitude of envelope [m]
% eta=real(a_0*exp(-((cg_0*(t-T/2)).^2/(2*sigma^2))+1i*omega_0*(t-T/2))*exp(1i*(phase_foc/180*pi)));

%%
sigma     = 1/delta/k_0; % group length of a gaussian packet is equal to 4 sigma
delta_x   = 1/sigma;         


%% double check
const_1   = omega_0^2*h_s/g;
k_0s      = fsolve(@(x) -const_1./x+tanh(x),const_1)/h_s; % wavenumber on the shallower side for given f_0 and h_S
k0shs     = k_0s*h_s;



%% make the shift of focus time due to the delay of wave generation
% this process makes the time the same with the experiments!
% T_start_WG = t_focus - (x_focus+7.5)/c_g0 ; 
% t_focus    = -T_start_WG + t_focus ; 

end

