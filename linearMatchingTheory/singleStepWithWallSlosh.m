clear
global zwuArr zwlArr
% close all
g = 9.81;
addpath c:/gits/wavemaker

%% plot and export
DO_EXPORT = true;
exportPath = './stepFigures/';
exportNamePrefix = 'grid_';
axisCommand = 'equal';
PLOT_X0_CONDITIONS = false;
PLOT_STREAMLINES = false;



ADD_OBSTACLE = 1;
RIGHT_WALL=1; xWall__L =.9;

L__lambda = .3;
T = 2;
w = 2*pi/T;
hL = 1.0;
hR = .8;
nEv = 400;
t = 0;
phaseDeg = 0;

lambda = 2*pi/findWaveNumbers(w,hR,0,0);
L = L__lambda*lambda;
xWall = xWall__L*L;

nx = 50;
dx = L/(nx-.5); % let x stop half a spacing befor x=0
xR = (.5:nx)*dx;
xL = -xR(nx:-1:1);
% xL = linspace(-1*L,0,50);
% xR = linspace(0,L,100);
z = linspace(-max([hL,hR]),0,100).';

dxPlot = hL*.003;

% obstacle
zwuArr = -[.2,.8]*hR;
zwlArr = -[.4,1]*hR; 

zwuArr = -[.6]*hR;
zwlArr = -[1]*hR; 

% grid
dz_o = 1/10*hR;
dz_w = .5*dz_o;
zwuArr = -(dz_o:2*dz_o:hR);
zwlArr = zwuArr-dz_w;


assert(all(zwlArr>=-hR));
assert(all(zwuArr-zwlArr>=0))

hphiI = 1.*exp( 1i* ( phaseDeg ) *pi/180  );


%% code
if ~ADD_OBSTACLE, [zwuArr,zwlArr]=deal(-hR); end


N = nEv+1;

kiR = findWaveNumbers(w,hR,0,nEv);
kiR = abs(real(kiR)) + 1i*(abs(imag(kiR)));
kiL = findWaveNumbers(w,hL,0,nEv);
kiL = abs(real(kiL)) + 1i*(abs(imag(kiL)));


kjR = [kiR;kiR].';
kjL = [kiL;kiL].';

wallTerm=0;
if RIGHT_WALL, wallTerm = exp(2i*kjR*xWall); end


A = nan(2*N);
b = nan(2*N,1);

d_i1 = zeros(N,1);d_i1(1) = 1;
[dL_j,dR_j] = deal(zeros(1,2*N));
dL_j(1:N) = 1;
dR_j(N+1:2*N) = 1;

dL_i  = [eye(N),  zeros(N)];
dCp_i = [zeros(N),eye(N), ];


% A:
A(1:N,:) = dL_i + dR_j.*(1-wallTerm).*kjR./kiL.*Gamma(kjR,hR,kiL,hL,'o')./Lambda(kiL,hL);
b(1:N) = hphiI.*d_i1;

% B:
A(N+1:2*N,:) = dR_j.*(1i*kiR.*(1-wallTerm).*Gamma(kjR,hR,kiR,hR,'w') + (1+wallTerm).*Gamma(kjR,hR,kiR,hR,'o')) ...
             - dL_j.*Gamma(kjL,hL,kiR,hR,'o');
b(N+1:2*N) = hphiI*Gamma(kjL(1),hL,kiR,hR,'o');

hphi = shiftdim(A\b,-2);


hphiL = hphi(:,:,1:N);
hphiR = hphi(:,:,N+1:2*N);
kiL = shiftdim(kiL,-2);
kiR = shiftdim(kiR,-2);



phiL = real( sum(  hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL.*xL-w*t)),3) +  hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(kiL(1).*xL-w*t)) );
%  phiR = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*(exp(1i*(kiR.*xR-w*t))+isWall*exp(1i*kiR.*(-xR+2*xWall)-1i*w*t)),3) );
phiR = real( sum(  hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(+kiR.*xR-w*t)).*(1+RIGHT_WALL*exp(-2i*kiR.*(xR-xWall))),3) );

if PLOT_STREAMLINES
    psiL = real( sum( -1i* hphiL.*sinh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL.*xL-w*t)),3) +  1i*hphiI.*sinh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(kiL(1).*xL-w*t)) );
    psiR = real( sum(  1i*hphiR.*sinh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(+kiR.*xR-w*t)).*(1-RIGHT_WALL*exp(-2i*kiR.*(xR-xWall))),3) );
    psiR(:,xR<0.025*hL)=nan;
    psiL(:,xL>-0.025*hL)=nan;
%     psiR(:,xR<0.05*lambda)=nan;
%     psiL(:,xL>-0.05*lambda)=nan;
end


x0L = -dxPlot;
x0R = dxPlot;

phiL_0 = real( sum( hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL*x0L-w*t)),3) + hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*(kiL(1)*x0L-w*t)) );
phiR_0 = real( sum( hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(kiR*x0R-w*t)).*(1+RIGHT_WALL*exp(-2i*kiR.*(x0R-xWall))),3) );

phiLx_0 = real( sum( -1i.*kiL.*hphiL.*cosh(kiL.*(hL+z))./cosh(kiL*hL).*exp(1i*(-kiL*x0L-w*t)),3) + 1i*kiL(1).* hphiI.*cosh(kiL(1).*(hL+z))./cosh(kiL(1)*hL).*exp(1i*( kiL(1)*x0L-w*t)) );
phiRx_0 = real( sum(  1i.*kiR.*hphiR.*cosh(kiR.*(hR+z))./cosh(kiR*hR).*exp(1i*(kiR.*x0R-w*t)).*(1-RIGHT_WALL*exp(-2i*kiR.*(x0R-xWall))),3) );
phiR_0(z<-hR)=nan;
phiRx_0(z<-hR)=nan;
%% plot


hf = figure('color','w','position',[1640 558 1129 420]);
if PLOT_X0_CONDITIONS,  subplot(1,3,1); end
[~,hc] = contourf([xL,xR],z,[phiL,phiR]);
hold on;
% if PLOT_STREAMLINES, contour([xL(xL<-0.05*hL),xR(xR>0.05*hL)],z,[psiL(:,xL<-0.05*hL),psiR(:,xR>0.05*hL)],'k'); end
if PLOT_STREAMLINES, contour([xL,xR],z,[psiL,psiR],'k'); end
ha = gca;
% plot([xL([1,end]),xR([1,end])],-[hL,hL,hR,hR],'k','linewidth',1.5)
patch([0,0,1,1]*xR(end),-[hL,hR,hR,hL],[1,1,1]*.5)
plot([0;0],[zwlArr;zwuArr],'k','linewidth',2)
if RIGHT_WALL
%     plot([1,1]*xWall,ha.YLim,'k-'); 
    patch([xWall*[1,1],ha.XLim(2)*[1,1]],-[hR,0,0,hR],[1,1,1]*.5)
end
axis(axisCommand)
% 
if PLOT_X0_CONDITIONS
    ha=subplot(1,3,2);
    plot(phiLx_0,z,'k',phiRx_0,z,'--r','linewidth',1);
    legend({'left','right'},'autoupdate','off')
    ylabel('z'); grid on;
    ha.XAxisLocation='origin';
    ha.YAxisLocation='origin';
    xlabel('\phi_x');
    hold on
    plot(ha.XLim',[-hR,zwuArr].*[1;1],'k--')
    %
    ha=subplot(1,3,3);
    plot(phiL_0,z,'k',phiR_0,z,'--r','linewidth',1);
    legend({'left','right'},'autoupdate','off')
    ylabel('z'); grid on;
    ha.XAxisLocation='origin';
    ha.YAxisLocation='origin';
    xlabel('\phi');
    hold on
    plot(ha.XLim',[-hR,zwuArr].*[1;1],'k--')
end

if DO_EXPORT
    if ~isfolder(exportPath), mkdir(exportPath); end
    fileName = sprintf('%sobst%d_wall%d_Lr%.2g_T%.2g_hL%.2g_hR%.2g_nEv%d',exportNamePrefix,ADD_OBSTACLE,RIGHT_WALL,L__lambda,T,hL,hR,nEv);
    fileName(fileName=='.')='p';
    
    copyfile('./singleStepWithWallSlosh.m',[exportPath,'/',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]) 
    export_fig(hf,[exportPath,'/',fileName],'-png','-pdf') 
end

function res = Lambda(k,h)
    res = .5*(h.*sech(k.*h).^2+tanh(k.*h)./k);
end
function res = Gamma(k1,h1,k2,h2,wallOrOpen)
    global zwuArr zwlArr
    
    switch wallOrOpen
        case 'o'
            zu = [0,zwlArr(1:end-1)];
            zl = zwuArr; 
            if zwlArr(end)>-min(h1,h2)
                zu(end+1) = zwlArr(end);
                zl(end+1) = -min(h1,h2);
            end
        case 'w'
            zu =  zwuArr;
            zl = zwlArr; 
            if h2>h1
                if zl(end)<=-h1
                    zl(end) = -h2;
                else
                    zu(end+1) = -h1;
                    zl(end+1) = -h2;
                end
            end
    end
    if all(zu==zl),res = 0; return;  end
    zu = shiftdim(zu,-1);    zl = shiftdim(zl,-1);
    iEq = k1==k2;
    if all(iEq,'all')
        res = sum(sech(k1*h1).*sech(k1*h2).*( .5*(zu-zl)*cosh(k1*(h1-h2)) + .25*(sinh(k1*(h1+h2+2*zu))-sinh(k1*(h1+h2+2*zl)))./k1 ),3);
        return
    end    
    res =  .5*sech(k1*h1).*sech(k2*h2).*sum(...
         ( sinh( k1.*(h1+zu) + k2.*(h2+zu) ) - sinh( k1.*(h1+zl) + k2.*(h2+zl) ) )./(k1+k2) ...
        +( sinh( k1.*(h1+zu) - k2.*(h2+zu) ) - sinh( k1.*(h1+zl) - k2.*(h2+zl) ) )./(k1-k2) ...
        ,3);
    if any(iEq,'all')
        k = k1+0*k2;
        k = k(iEq);
        res(iEq) = sech(k*h1).*sech(k*h2).*sum( .5*(zu-zl).*cosh(k*(h1-h2)) + .25*(sinh(k.*(h1+h2+2*zu))-sinh(k.*(h1+h2+2*zl)))./k ,3);
    end
end

