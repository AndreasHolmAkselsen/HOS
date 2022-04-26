clear
% close all
g = 9.81;
addpath c:/gits/wavemaker

% T = 4.1733;
T = 2.2987;
w = 2*pi/T;
hL = 1.0;
hR = .5;
nEv = 500;
t = 0;

hphiI = exp(1i*deg2rad(90));


N = nEv+1;

kiR = findWaveNumbers(w,hR,0,nEv);
kiR = abs(real(kiR)) + 1i*(abs(imag(kiR)));
kiL = findWaveNumbers(w,hL,0,nEv);
kiL = abs(real(kiL)) + 1i*(abs(imag(kiL)));


L = .25*pi/kiL(1);
xL = linspace(-1*L,0,50);
xR = linspace(0,L,100); xR(1)= [];
z = linspace(-max([hL,hR]),0,90).';

kjR = [kiR;kiR].';
kjL = [kiL;kiL].';

A = nan(2*N);
b = nan(2*N,1);

d_i1 = zeros(N,1);d_i1(1) = 1;
[dL_j,dRp_j] = deal(zeros(1,2*N));
dL_j(1:N) = 1;
dRp_j(N+1:2*N) = 1;




% A:
% dL_i  = [eye(N),  zeros(N)];
% A(1:N,:) = dL_i + dRp_j.* kjR./kiL.*Gamma(kiL,hL,kjR,hR)./Lambda(kiL,hL);

A(1:N,1:N) = eye(N);
A(1:N,N+1:2*N) = kiR.'./kiL.*Gamma(kiL,hL,kiR.',hR)./Lambda(kiL,hL);

b(1:N) = hphiI.*d_i1;


% B:
% dRp_i = [zeros(N),eye(N), ];
% A(N+1:2*N,:) = dRp_i - dL_j.*Gamma(kjL,hL,kiR,hR)./Lambda(kiR,hR);
A(N+1:2*N,N+1:2*N) = eye(N);
A(N+1:2*N,1:N) = -Gamma(kiL.',hL,kiR,hR)./Lambda(kiR,hR);

b(N+1:2*N) = hphiI*Gamma(kjL(1),hL,kiR,hR)./Lambda(kiR,hR);

% A./b


hphi = shiftdim(A\b,-2);

hphiL = hphi(:,:,1:N);
hphiR = hphi(:,:,N+1:2*N);
kiL = shiftdim(kiL,-2);
kiR = shiftdim(kiR,-2);
phiL = real( sum(  hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL.*xL-w*t)),3) +  hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(kiL(1).*xL-w*t)) );
phiR = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(+kiR.*xR-w*t)),3) );

phiL_0 = real( sum( hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3) + hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiR_0 = real( sum( hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3) );

phiLx_0 = real( sum( -1i.*kiL.*hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3) + 1i*kiL(1).* hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiRx_0 = real( sum( 1i*kiR.* hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3) );

phiR_0(z<-hR)=nan;
phiRx_0(z<-hR)=nan;


fprintf('Reflection coefficient: %g\n', abs(hphiL(1))/abs(hphiI) )

%% plot


hf = figure('color','w','position',[1640 558 1129 240]);
subplot(1,3,1);
contourf([xL,xR],z,[phiL,phiR]);
hold on
patch([0,0,L,L],-[hL,hR,hR,hL],.5*[1,1,1],'lineStyle','none');%,'FaceAlpha',.5
plot([-L,0,0,L],-[hL,hL,hR,hR],'k','linewidth',1.)
title('\phi')
xlabel('x');ylabel('y')
axis equal
% 
ha=subplot(1,3,2);
plot(phiLx_0,z,'k',phiRx_0,z,'--r','linewidth',1);
xlabel('\phi_x'); ylabel('y'); grid on;
% ha.XAxisLocation='origin';
% 
ha=subplot(1,3,3);
plot(phiL_0,z,'k',phiR_0,z,'--r','linewidth',1);
xlabel('\phi'); ylabel('y'); grid on;
% ha.XAxisLocation='origin';


return
export_fig('contourPlot','-pdf','-png','-m2')

function res = Lambda(k,h)
    res = .5*(h.*sech(k.*h).^2+tanh(k.*h)./k);
end
function res = Gamma(k1,h1,k2,h2)
    res =  (k1.*tanh(k1.*h1)-k2.*tanh(k2.*h2)-k1.*sinh(k1*(h1-h2)).*sech(k1.*h1).*sech(k2.*h2))./(k1.^2-k2.^2);
    ii = k1==k2;
    if any(ii,'all')
        k = k1+0*k2;
        k = k(ii);
        res(ii) = .5*h2 + .5*(1-k.*h2.*tanh(k*h1))./k.*tanh(k*h2);
    end
%     res = res./(.5*(h2.*sech(k2.*h2).^2+tanh(k2.*h2)./k2));
end


% % eq tests:
% hphiR_beta = -hphiL + ((1:N)'==1)*hphiI;
% 
% k = [-k;k];
% hphiL_tot = [hphiI;zeros(nEv,1);hphiL]; % [kk; +ikk,-kk;-ikk];
% hphiR_tot = [hphiR_beta;zeros(nEv+1,1)];
% 
% mk = [(1:N)+N,1:N]';
% eq_R = k.*hphiL_tot + k(mk).* hphiL_tot(mk) - ( k.* hphiR_tot + k(mk) .*hphiR_tot(mk)  );
% assert(max(abs(eq_R./(k(1)*hphiI))) < 1e-9 )
% 
% 
% l = k.';
% Lambda = .5*(h.*sech(k.*h).^2+tanh(k.*h)./k);
% Gamma = sech(k*h).*sech(l*h) .*( .5*(sinh((k+l).*(h+zu)) - sinh((k+l).*(h+zl)))./(k+l) + .5*(sinh((k-l).*(h+zu)) - sinh((k-l).*(h+zl)))./(k-l)  );
% Gamma(k==l) = sech(k*h).^2 .*( .5*(zu-zl) + .25./k.*(sinh(2*k.*(h+zu))-sinh(2*k.*(h+zl))) );
% Gamma(k==-l) = Gamma(k==l);
% a_tot = Gamma./Lambda;
% 
% eq_A = hphiL_tot + hphiL_tot(mk) - sum(a_tot.*hphiR_tot.',2);
% assert(max(abs(eq_A./hphiI)) < 1e-9)

% % eq_B = hphiR_tot - hphiR_tot(mk) - sum(a_tot.*hphiL_tot.',2)
