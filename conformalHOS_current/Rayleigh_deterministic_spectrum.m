clear


N_ensamble= 200;
L = 100;
nx = 5000;
Tp = 2;
Hs = .1;
gamma = 3;
T = 1000;


df = 1/T;
nf = floor(nx/2)-1;
f = (0:nf-1)'*df;
S = level1.wave.computeWaveSpectralDensity_JONSWAP(f,Hs,Tp,gamma);

R = @(n,m) (randn(n,m) + 1i*randn(n,m))/sqrt(2);
% abs(R) is Rayleigh-distributed and angle(R) is uniformly distributed.
% Rrand is scaled such that the second moment E[Rrand^2] = 1.
% Per definition, S=E[lim_T 1/T*|hx|^2] -> lim_m lim_T a^2/T *|R_m|^2
% : hx = a*Rrand with a = sqrt(S*T)
%        test mean(abs(R(n,1).^2)) -> 1.0

% deterministic amplitudes:  R = @(n,m) exp(2i*pi*rand(n,m));

hat_x_plus = sqrt(S*T).*R(nf,N_ensamble);  % integral definition of modes.       
S_rand = mean(abs(hat_x_plus_rand).^2 /T ,2); % arbitrary T is added and here devided away
figure, plot(f,S,f,S_rand,'.')

% for discrete modes (sum instead of integral) T is replaced by df above.
% See Næss&Moan eq. (6.35), example 6.3.1 and section 8.2.1.