clear
close all
g = 9.81;
addpath c:/gits/wavemaker


%% plot and export
DO_EXPORT = false;
exportPath = './doubleStepFigures/';
exportNamePrefix = '';
axisCommand = 'equal';
PLOT_X0_CONDITIONS = false;

L__lambda = .3;
T = 1.5;
w = 2*pi/T;
hL = .9;
hC = .5;
hR = .8;

lambda = 2*pi/findWaveNumbers(w,hR,0,0);
L = L__lambda*lambda;


nEv = 400;
t =  0*5643.353;
epsx = .001*hL;
xL = linspace(-1,0-epsx,100);
xC = linspace(0+epsx,L-epsx,200); xC(1)= [];
xR = linspace(L+epsx,L+1,100); xR(1)= [];
z = linspace(-max([hL,hC,hR]),0,200).';

hphiI = 1.*exp( 1i* ( 0 ) *pi/180  );

N = nEv+1;

kiR = findWaveNumbers(w,hR,0,nEv);
kiR = abs(real(kiR)) + 1i*(abs(imag(kiR)));
kiC = findWaveNumbers(w,hC,0,nEv);
kiC = abs(real(kiC)) + 1i*(abs(imag(kiC)));
kiL = findWaveNumbers(w,hL,0,nEv);
kiL = abs(real(kiL)) + 1i*(abs(imag(kiL)));


kjR = [kiR;kiR;kiR;kiR].';
kjC = [kiC;kiC;kiC;kiC].';
kjL = [kiL;kiL;kiL;kiL].';

A = nan(4*N);
b = nan(4*N,1);


d_i1 = zeros(N,1);d_i1(1) = 1;
[dL_j,dCp_j,dCm_j,dR_j] = deal(zeros(1,4*N));
dL_j(1:N) = 1;
dCp_j(N+1:2*N) = 1;
dCm_j(2*N+1:3*N) = 1;
dR_j(3*N+1:4*N) = 1;

dL_i  = [eye(N),  zeros(N),zeros(N),zeros(N)];
dCp_i = [zeros(N),eye(N),  zeros(N),zeros(N)];
dCm_i = [zeros(N),zeros(N),eye(N),  zeros(N)];
dR_i  = [zeros(N),zeros(N),zeros(N),eye(N)  ];


% A:
A(1:N,:) = dL_i + (dCp_j-dCm_j.*exp(1i*kjC*L)).* kjC./kiL.*Gamma(kiL,hL,kjC,hC)./Lambda(kiL,hL);
b(1:N) = hphiI.*d_i1;

% B:
A(N+1:2*N,:) = dCp_i+dCm_i.*exp(1i*kiC*L) - dL_j.*Gamma(kjL,hL,kiC,hC)./Lambda(kiC,hC);
b(N+1:2*N) =  hphiI*Gamma(kiL(1),hL,kiC,hC)./Lambda(kiC,hC);

% C:
A(2*N+1:3*N,:) = dR_i - (dCp_j.*exp(1i*kjC*L)-dCm_j).*kjC./kiR.*Gamma(kiR,hR,kjC,hC)./Lambda(kiR,hR);
b(2*N+1:3*N) = 0;

% D:
A(3*N+1:4*N,:) = (dCp_i.*exp(1i*kiC*L)+dCm_i) - dR_j.*Gamma(kjR,hR,kiC,hC)./Lambda(kiC,hC);
b(3*N+1:4*N) = 0;

hphi = shiftdim(A\b,-2);

hphiL = hphi(:,:,1:N);
hphiCp = hphi(:,:,N+1:2*N);
hphiCm = hphi(:,:,2*N+1:3*N);
hphiR = hphi(:,:,3*N+1:4*N);
kiL = shiftdim(kiL,-2);
kiC = shiftdim(kiC,-2);
kiR = shiftdim(kiR,-2);
phiL = real( sum(  hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL.*xL-w*t)),3) +  hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(kiL(1).*xL-w*t)) );
phiC = real( sum(  hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(+kiC.*xC-w*t)),3) + sum(  hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-kiC.*(xC-L)-w*t)),3) );
phiR = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(+kiR.*(xR-L)-w*t)),3));

phiLx_0 = real( sum( -1i.*kiL.*hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3) + 1i*kiL(1).* hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiCx_0 = real( sum( 1i*kiC.* hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3) + sum( -1i*kiC.* hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-kiC.*(-L)-w*t)),3) );
phiCx_L = real( sum( 1i*kiC.*hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(+kiC.*L-w*t)),3) + sum(  -1i*kiC.*hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3) );
phiRx_L = real( sum( 1i*kiR.* hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3));

phiL_0 = real( sum( hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-w*t)),3) +  hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(-w*t)) );
phiC_0 = real( sum( hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3) + sum(  hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-kiC.*(-L)-w*t)),3) );
phiC_L = real( sum( hphiCp.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(+kiC.*L-w*t)),3) + sum( hphiCm.*cosh(kiC.*(hC+z))./cosh(kiC*hC).*exp(1i*(-w*t)),3) );
phiR_L = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(-w*t)),3));

phiLx_0(z<-hL)=nan;
phiCx_0(z<-hC)=nan;
phiCx_L(z<-hC)=nan;
phiRx_L(z<-hR)=nan;
phiL_0(z<-hL)=nan;
phiC_0(z<-hC)=nan;
phiC_L(z<-hC)=nan;
phiR_L(z<-hR)=nan;
%% plot




hf = figure('color','w','position',[1640 558 1129 420]);
if PLOT_X0_CONDITIONS,  subplot(1,3,1); end
contourf([xL,xC,xR],z,[phiL,phiC,phiR]);
hold on
% plot([xL([1,end]),xC([1,end]),xR([1,end])],-[hL,hL,hC,hC,hR,hR],'k','linewidth',1.5)
patch([xL(1),xL(1),0,0,L,L,xR(end),xR(end)],-[1,hL,hL,hC,hC,hR,hR,1],[1,1,1]*.5)

if PLOT_X0_CONDITIONS
    ha=subplot(1,3,2);
    plot(phiLx_0,z,'k',phiCx_0,z,'--r',phiCx_L,z,'g',phiRx_L,z,'--b','linewidth',1);
    ylabel('z'); grid on;
    ha.XAxisLocation='origin';
    xlabel('\phi_x');
    
    ha=subplot(1,3,3);
    plot(phiL_0,z,'k',phiC_0,z,'--r',phiC_L,z,'g',phiR_L,z,'--b','linewidth',1);
    ylabel('z'); grid on;
    ha.XAxisLocation='origin';
    xlabel('\phi');
end



if DO_EXPORT
    if ~isfolder(exportPath), mkdir(exportPath); end
    fileName = sprintf('%s_Lr%.2g_T%.2g_hL%.2g_hC%.2g_hR%.2g_nEv%d',exportNamePrefix,L__lambda,T,hL,hC,hR,nEv);
    fileName(fileName=='.')='p';
    
    copyfile('./doubleStep.m',[exportPath,'/',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]) 
    export_fig(hf,[exportPath,'/',fileName],'-png','-pdf') 
end


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


