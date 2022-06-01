function [h,phiS] = initSSGW(k0,H,ka,N_SSGW,NWaves,x,g)

[zIC,dwdz,PP] = SSGW(k0.*H,ka,N_SSGW);

if isinf(PP(1)), L_scale = 1/k0; else, L_scale = H; end
out.c_e = PP(4)*sqrt(g*L_scale); % phase velocity observed from where the meam velocity at the bed is zero
out.c_s = PP(5)*sqrt(g*L_scale); % mean flow velocity (phase velocity in frame without mean flow)
out.k = PP(2)/L_scale;
zIC = zIC*L_scale;

%         % to move wave to centre (optional)
%         z = [ z(N_SSGW+1:end)-lambda/2 ; z(1:N_SSGW)+lambda/2 ];
%         dwdz = [ dwdz(N_SSGW+1:end); dwdz(1:N_SSGW) ];

% duplicate across domain.
lambda = 2*pi/k0;
zIC = reshape(repmat(zIC,1,NWaves)+lambda*(0:NWaves-1),[],1) + x(1);
dwdz = repmat(dwdz,NWaves,1);
dwdz = dwdz*sqrt(g*L_scale);

n = 2*N_SSGW*NWaves;
% z_m = .5*(zIC(1:n-1)+zIC(2:n));
dwdz0_m = .5*(dwdz(1:n-1)+dwdz(2:n))+out.c_e;
wIC = [0;cumsum( dwdz0_m.*diff(zIC))];
wIC = wIC-mean(wIC);

% if z(1)<2*eps&&z(1)>-2*eps, z(1)=1i*imag(z(1));end
Lx = -2*x(1);
zIC = [zIC(end)-Lx;zIC;zIC(1)+Lx]; % extend with ghost nodes
wIC = [wIC(end);wIC;wIC(end)];
h = interp1(real(zIC),imag(zIC),x,'linear',nan);
phiS = interp1(real(zIC),real(wIC),x,'linear',nan);
fft_h = fftshift(fft(h));
if sum(abs(fft_h(1:floor(end/4))))>.01*sum(abs(fft_h))
    warning('Initial conition may not have been found. Verify that solution exists.')
end
%         phiS00 = ka/k0.*g/omega*sin(xk0-phaseAng);
%         h00 = ka/k0*(cos(xk0-phaseAng));
%         % phi = ka/k0.*g/omega*sin(xk0-phaseAng)*cosh(k*(H+z))/cosh(k*H);
%         u0 = ka/k0.*g/omega*k0*cos(xk0-phaseAng);
%         v0 =  ka/k0.*g/omega*k0*sin(xk0-phaseAng)*tanh(k0*H_IC);
%         figure('color','w')
%         subplot(311), plot(x,h0,'-',x,h00,'--');ylabel('\eta'); grid on; legend('SSGW','linear')
%         subplot(312), plot(x,phiS0,'-',x,phiS00,'--');ylabel('\phi^S'); grid on; legend('SSGW','linear')
%         subplot(313), plot(x,u0,'-r',x,v0,'-b',real(zIC),real(dwdz)+out.c_e,'--r',real(zIC),-imag(dwdz),'--b');ylabel('velocity'); grid on
