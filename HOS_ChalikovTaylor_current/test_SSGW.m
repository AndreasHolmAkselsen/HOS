close all
clear
DO_EXPORT = false;

g = 9.81;

h = 0.5;
lambda = 10;
H = .025*lambda;



N = 500;
k = 2*pi/lambda;
kh = k*h;
[z,U,PP] = SSGW(kh,k*H/2,N);
out.hk = PP(1)*PP(2);
out.c_e = PP(4)*sqrt(g*h); % phase velocity observed from where the meam velocity at the bed is zero
out.c_s = PP(5)*sqrt(g*h); % mean flow velocity (phase velocity in frame without mean flow)

z = [ z(N+1:end)-2*pi/kh ; z(1:N) ];
U = [ U(N+1:end); U(1:N) ];
n = length(z);
z_m = .5*(z(1:n-1)+z(2:n));
U0_m = .5*(U(1:n-1)+U(2:n))+out.c_e;
w = [0;cumsum( conj(U0_m).*diff(z))];
w = w-mean(w);
%   In deep water:   rho = g = k = 1.
%   In finite depth: rho = g = d = 1.
%       PP(1)=depth, PP(2)=wavenumber, PP(3)=waveheight, 
%       PP(4)=celerity c_e, PP(5)=celerity c_s, PP(6)=Bernoulli constant, 
%       PP(7)=crest height, PP(8)=trough height, PP(9)=impulse, 
%       PP(10)=potential energy, PP(11)=kinetic energy, PP(12)=radiation stress,
%       PP(13)=momentum flux, PP(14)=energy flux, PP(15)=group velocity.


hf = figure('color','w');%,'Position',[1.3523   -0.0617    0.8500    0.5000]*1e3);
subplot(2,1,1); plot(z);ylabel('\eta');grid on;
subplot(2,1,2); plot(real(z),real(w));ylabel('\phi^S');grid on;


if DO_EXPORT
    savefig(hf,['.\results\SSGW\',TLabel]);
    export_fig(hf,['.\results\SSGW\',TLabel],'-png','-pdf')
end