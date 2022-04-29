clear
close all
g = 9.81;
addpath c:/gits/wavemaker

L = 1.3;
T = 2;
w = 2*pi/T;
hC = .7;
hR = .9;

nEv = 200;
t =  0*5643.353;
xC = linspace(0,L,100); xC(1)= [];
xR = linspace(L,L+1,50); xR(1)= [];
z = linspace(-max([hC,hR]),0,90).';

hphiI = 1;


N = nEv+1;

kiR = findWaveNumbers(w,hR,0,nEv);
kiR = abs(real(kiR)) + 1i*(abs(imag(kiR)));
kiC = findWaveNumbers(w,hC,0,nEv);
kiC = abs(real(kiC)) + 1i*(abs(imag(kiC)));


kjR = [kiR;kiR].';
kjC = [kiC;kiC].';

A = nan(2*N);
b = nan(2*N,1);


d_i1 = zeros(N,1);d_i1(1) = 1;
[dL_j,dCp_j,dCm_j,dR_j] = deal(zeros(1,2*N));
dCm_j(1:N) = 1;
dR_j(N+1:2*N) = 1;

dCm_i = [eye(N),  zeros(N)];
dR_i  = [zeros(N),eye(N)  ];

fixSwitch = 1;

% C:
A(1:N,:) = dR_i - (-dCm_j).*kjC./kiR.*Gamma(kiR,hR,kjC,hC)./Lambda(kiR,hR);
b(1:N) =  hphiI.*d_i1;

% D:
A(N+1:2*N,:) = (+dCm_i) - dR_j.*Gamma(kjR,hR,kiC,hC)./Lambda(kiC,hC);
b(N+1:2*N) = -fixSwitch*hphiI*d_i1;

hphi = shiftdim(A\b,-2);

hphiCm = hphi(:,:,1:N);
hphiR = hphi(:,:,N+1:2*N);
kiC = shiftdim(kiC,-2);
kiR = shiftdim(kiR,-2);
phiC = real( sum(  hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-kiC.*(xC-L)-w*t)),3)   + fixSwitch* hphiI.*cosh(kiC(1).*(hC+z))./cosh(kiC(1)*hC).*exp(1i*(kiC(1).*(xC-L)-w*t))  );
phiR = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(+kiR.*(xR-L)-w*t)),3));

phiCx_L = real( sum(  -1i*kiC.*hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3)   +fixSwitch* 1i*kiC(1)*hphiI.*cosh(kiC(1).*(hC+z))./cosh(kiC(1)*hC).*exp(1i*(-w*t)) );
phiRx_L = real( sum( 1i*kiR.* hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3));

phiC_L = real(sum( hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3)+fixSwitch* hphiI.*cosh(kiC(1).*(hC+z))./cosh(kiC(1)*hC).*exp(1i*(-w*t)) );
phiR_L = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3));

phiCx_L(z<-hC)=nan;
phiRx_L(z<-hR)=nan;
phiC_L(z<-hC)=nan;
phiR_L(z<-hR)=nan;
%% plot




hf = figure('color','w','position',[1640 558 1129 420]);
subplot(1,3,1);
contourf([xC,xR],z,[phiC,phiR]);
hold on
plot([xC([1,end]),xR([1,end])],-[hC,hC,hR,hR],'k','linewidth',1.5)


% 
ha=subplot(1,3,2);
plot(phiCx_L,z,'g',phiRx_L,z,'--b','linewidth',1);
ylabel('z'); grid on;
ha.XAxisLocation='origin';
xlabel('\phi_x');

ha=subplot(1,3,3);
plot(phiC_L,z,'g',phiR_L,z,'--b','linewidth',1);
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
% eq_C = k.*hphiL_tot + k(mk).* hphiL_tot(mk) - ( k.* hphiR_tot + k(mk) .*hphiR_tot(mk)  );
% assert(max(abs(eq_C./(k(1)*hphiI))) < 1e-9 )
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
