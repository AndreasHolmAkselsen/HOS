clear
close all
g = 9.81;
addpath c:/gits/wavemaker


L = 4;
T = 4.1733;
w = 2*pi/T;
hR = 1.0;
hL = .5;
nEv = 400;

t =  0*5643.353;
xL = linspace(0,L,100); xL(1)= [];
xR = linspace(L,L+1,50); xR(1)= [];
z = linspace(-max([hL,hR]),0,90).';

hphiI = 1;


N = nEv+1;

kiR = findWaveNumbers(w,hR,0,nEv);
kiR = abs(real(kiR)) + 1i*(abs(imag(kiR)));
kiL = findWaveNumbers(w,hL,0,nEv);
kiL = abs(real(kiL)) + 1i*(abs(imag(kiL)));


kjR = [kiR;kiR].';
kjL = [kiL;kiL].';

A = nan(2*N);
b = nan(2*N,1);


d_i1 = zeros(N,1);d_i1(1) = 1;
[dL_j,dLp_j,dLm_j,dR_j] = deal(zeros(1,2*N));
dLm_j(1:N) = 1;
dR_j(N+1:2*N) = 1;

dLm_i = [eye(N),  zeros(N)];
dR_i  = [zeros(N),eye(N)  ];

fixSwitch = 1;

% L:
A(1:N,:) = dR_i - (-dLm_j).*kjL./kiR.*Gamma(kiR,hR,kjL,hL)./Lambda(kiR,hR);
b(1:N) =  hphiI.*d_i1;

% D:
A(N+1:2*N,:) = (+dLm_i) - dR_j.*Gamma(kjR,hR,kiL,hL)./Lambda(kiL,hL);
b(N+1:2*N) = -fixSwitch*hphiI*d_i1;

hphi = shiftdim(A\b,-2);

hphiL = hphi(:,:,1:N);
hphiR = hphi(:,:,N+1:2*N);
kiL = shiftdim(kiL,-2);
kiR = shiftdim(kiR,-2);
phiL = real( sum(  hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL.*(xL-L)-w*t)),3)   + fixSwitch* hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(kiL(1).*(xL-L)-w*t))  );
phiR = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(+kiR.*(xR-L)-w*t)),3));

phiLx_L = real( sum(  -1i*kiL.*hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3)   +fixSwitch* 1i*kiL(1)*hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiRx_L = real( sum( 1i*kiR.* hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3));

phiL_L = real(sum( hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3)+fixSwitch* hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiR_L = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3));

phiLx_L(z<-hL)=nan;
phiRx_L(z<-hR)=nan;
phiL_L(z<-hL)=nan;
phiR_L(z<-hR)=nan;
%% plot


fprintf('Reflection coefficient: %g\n', abs(hphiL(1))/abs(hphiI) )


hf = figure('color','w','position',[1640 558 1129 420]);
subplot(1,3,1);
contourf([xL,xR],z,[phiL,phiR]);
hold on
plot([xL([1,end]),xR([1,end])],-[hL,hL,hR,hR],'k','linewidth',1.5)


% 
ha=subplot(1,3,2);
plot(phiLx_L,z,'g',phiRx_L,z,'--b','linewidth',1);
ylabel('z'); grid on;
ha.XAxisLocation='origin';
xlabel('\phi_x');

ha=subplot(1,3,3);
plot(phiL_L,z,'g',phiR_L,z,'--b','linewidth',1);
ylabel('z'); grid on;
ha.XAxisLocation='origin';
xlabel('\phi');

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
% eq_L = k.*hphiL_tot + k(mk).* hphiL_tot(mk) - ( k.* hphiR_tot + k(mk) .*hphiR_tot(mk)  );
% assert(max(abs(eq_L./(k(1)*hphiI))) < 1e-9 )
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
