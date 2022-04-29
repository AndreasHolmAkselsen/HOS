clear
% close all
g = 9.81;
addpath c:/gits/wavemaker

L = 1;
T = 2;
w = 2*pi/T;
hL = 1.0;
hC = .6;
nEv = 200;
t = 543.353;
xL = linspace(-1*L,0,50);
xC = linspace(0,L,100); xC(1)= [];
z = linspace(-max([hL,hC]),0,90).';

hphiI = 1+.76i;


N = nEv+1;

kiC = findWaveNumbers(w,hC,0,nEv);
kiC = abs(real(kiC)) + 1i*(abs(imag(kiC)));
kiL = findWaveNumbers(w,hL,0,nEv);
kiL = abs(real(kiL)) + 1i*(abs(imag(kiL)));


kjC = [kiC;kiC].';
kjL = [kiL;kiL].';

A = nan(2*N);
b = nan(2*N,1);

d_i1 = zeros(N,1);d_i1(1) = 1;
[dL_j,dCp_j] = deal(zeros(1,2*N));
dL_j(1:N) = 1;
dCp_j(N+1:2*N) = 1;




% A:
% dL_i  = [eye(N),  zeros(N)];
% A(1:N,:) = dL_i + dCp_j.* kjC./kiL.*Gamma(kiL,hL,kjC,hC)./Lambda(kiL,hL);

A(1:N,1:N) = eye(N);
A(1:N,N+1:2*N) = kiC.'./kiL.*Gamma(kiL,hL,kiC.',hC)./Lambda(kiL,hL);

b(1:N) = hphiI.*d_i1;


% B:
% dCp_i = [zeros(N),eye(N), ];
% A(N+1:2*N,:) = dCp_i - dL_j.*Gamma(kjL,hL,kiC,hC)./Lambda(kiC,hC);
A(N+1:2*N,N+1:2*N) = eye(N);
A(N+1:2*N,1:N) = -Gamma(kiL.',hL,kiC,hC)./Lambda(kiC,hC);

b(N+1:2*N) = hphiI*Gamma(kjL(1),hL,kiC,hC)./Lambda(kiC,hC);

% A./b


hphi = shiftdim(A\b,-2);

hphiL = hphi(:,:,1:N);
hphiCp = hphi(:,:,N+1:2*N);
kiL = shiftdim(kiL,-2);
kiC = shiftdim(kiC,-2);
phiL = real( sum(  hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL.*xL-w*t)),3) +  hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(kiL(1).*xL-w*t)) );
phiC = real( sum(  hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(+kiC.*xC-w*t)),3) );

phiL_0 = real( sum( hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3) + hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiC_0 = real( sum( hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3) );

phiLx_0 = real( sum( -1i.*kiL.*hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3) + 1i*kiL(1).* hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiCx_0 = real( sum( 1i*kiC.* hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3) );

phiC_0(z<-hC)=nan;
phiCx_0(z<-hC)=nan;
%% plot


hf = figure('color','w','position',[1640 558 1129 420]);
subplot(1,3,1);
contourf([xL,xC],z,[phiL,phiC]);
hold on
plot([xL([1,end]),xC([1,end])],-[hL,hL,hC,hC],'k','linewidth',1.5)
% 
ha=subplot(1,3,2);
plot(phiLx_0,z,'k',phiCx_0,z,'--r','linewidth',1);
ylabel('z'); grid on;
ha.XAxisLocation='origin';
xlabel('\phi_x');
% 
ha=subplot(1,3,3);
plot(phiL_0,z,'k',phiC_0,z,'--r','linewidth',1);
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
