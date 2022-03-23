clear
global M x k_cut nonLinRamp  surfaceMethod timeReached t_end H dW DO_PADDING kx
timeReached = 0;

%% input
% surfaceMethod = 'decayingConformal'; % Chalikov method
surfaceMethod = 'Taylor';  % normal HOS
DO_PADDING = true;

% Resolution
nx__wave = 2^10/10;
M = 5; % solution order for Taylor method
relTolODE = 1e-4;% 1e-8;
N_SSGW = 2000; % number of modes in SSGW solution

% Plot & export options
DO_EXPORT = 1;
EXPORT_MAT = 1;
PLOT_CURRENT = false;
exportPrefix = '';
exportPath = './figures/';
i_detailedPlot = 5; %plot contour plots of frame i. Leave empty to skip

% Wave specification
NWaves = 10;
lambda = 10;
ka = .4; % linear wave steepness

h = 100;%2*lambda; % water depth. 

% some computations...
L = NWaves*lambda;
g = 9.81;
k0 = 2*pi/lambda;
% omega = (1+.5*ka^2)*sqrt(g*k0*tanh(k0*H));
omega = sqrt(g*k0*tanh(k0*h));
T = 2*pi/omega;
c_p = 2*pi/T/k0;
nx = nx__wave*NWaves;


% Simulation/plotting time
NT_dt = 2.5;
dt = NT_dt*T;
t_end = 9*dt;

    
% stability
% DO_LIN_WAVE_INIT = false;
% Tramp = 0*1*T;

DO_LIN_WAVE_INIT = true;
Tramp = 10*T;

nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
% k_cutTaylor = (M+5)*k0;
k_cutTaylor = 50*(2*pi/L);
k_cut_conformal =  (nx*pi/L)/2;



%% Simulation
dx = L/nx;
x = (0:nx-1)'*dx;% - L/2;
t0 = 0;
initialStepODE = 1e-3*T;
xk0 = k0.*x;
phaseAng = 0*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);

%% Specify background current
zeta_j = []; F_j = [];
nMirror = 3; % number of times the domain is repeated in x.

% single vortex
zeta_j = [.5-.075i  ]*L;% object centre
F_j    = 0*[ -.2i  ];% object strength


F_j = shiftdim(F_j,-1); zeta_j = shiftdim(zeta_j,-1);% ID_j = shiftdim(ID_j,-1);
zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);


% % vortex/source/sink
A_j = .5*F_j.*c_p.*abs(imag(zeta_j));
W  = @(zz) sum(A_j.*log(zz-zeta_j) + conj(A_j.*log(conj(zz)-zeta_j)),3:4);
dW = @(zz) sum(A_j./(zz-zeta_j) + conj(A_j./(conj(zz)-zeta_j)),3:4);

% doublet
% A_j = .5*F_j.*c_p.*abs(imag(zeta_j)).^2;
% f = @(zz) sum(-A_j./(zz-zeta_j) - conj(A_j)./(zz-conj(zeta_j)),3:4);
% dW = @(zz) sum(A_j./(zz-zeta_j).^2 + conj(A_j)./(zz-conj(zeta_j)).^2,3:4);

if all(F_j==0), dW=@(zz)0; end



if PLOT_CURRENT && ~isempty(F_j)
    z = 0:dx:.3*L; z = [-z(end:-1:2),z];
    hf_c = figure('color','w'); ha = gca;
    % plot velocity intensity |U|
    absU = abs(dW(x+1i*z));
    ULim = 2.5*max(abs(dW(x)));
    hIm = imagesc(x',z',absU',[0,ULim]); 
    hIm.AlphaData = (absU<ULim)';
    % plot streamlines
    hold on, axis equal xy
    contour(x',z',imag(W(x+1i*z))',20,'k');
    plot(ha.XLim,[0,0],'k')
    ylim([z(1),0])
    
    % add quiver plot
    zz_ip = linspace(x(1),x(end),10) + 1i*linspace(z(1),z(end),10)';
    dW_ip = dW(zz_ip);
    dW_ip(abs(dW_ip)>ULim) = nan;
    quiver(real(zz_ip),imag(zz_ip),real(dW_ip),-imag(dW_ip),'r');
    
    ha.Visible='off';
%     fprintf('c_p = %.3gm/s, max |U(0)| = %.3gm/s, fraction: %.3g\n',c_p,max(abs(dW(x))),max(abs(dW(x)))/c_p)
    drawnow
    
    if DO_EXPORT
        fileName = ['curr_',exportPrefix];
        savefig(hf_c,[exportPath,'/',fileName]);
        export_fig(hf_c,[exportPath,'/',fileName],'-pdf');
    end
    
end


%% init with SSGW:

if DO_LIN_WAVE_INIT
    eta = ka/k0*(cos(xk0-phaseAng));
    phiS = ka/k0.*g/omega*sin(xk0-phaseAng);
else
    
    [z,dwdz,PP] = SSGW(k0*h,ka,N_SSGW);
    

    
    
    % z = [ z(N_SSGW+1:end)-2*pi/(k0*h) ; z(1:N_SSGW) ];
    % z = z*h;
    % dwdz = [ dwdz(N_SSGW+1:end); dwdz(1:N_SSGW) ];
    % z = reshape(repmat(z,1,NWaves) + L*(floor(-NWaves/2+1):floor(NWaves/2)),[],1);
    
    if isinf(PP(1)), L_scale = 1/k0; else, L_scale = h; end
    
    
    out.c_e = PP(4)*sqrt(g*L_scale); % phase velocity observed from where the meam velocity at the bed is zero
    out.c_s = PP(5)*sqrt(g*L_scale); % mean flow velocity (phase velocity in frame without mean flow)
    out.k = PP(2)/L_scale;
    z = z*L_scale;
    
    z = reshape(repmat(z,1,NWaves)+lambda*(0:NWaves-1),[],1);
    dwdz = repmat(dwdz,NWaves,1);
    dwdz = dwdz*sqrt(g*L_scale);
    
    n = 2*N_SSGW*NWaves;
    z_m = .5*(z(1:n-1)+z(2:n));
    dwdz0_m = .5*(dwdz(1:n-1)+dwdz(2:n))+out.c_e;
    w = [0;cumsum( dwdz0_m.*diff(z))];
    w = w-mean(w);
    
    % if z(1)<2*eps&&z(1)>-2*eps, z(1)=1i*imag(z(1));end
    z_ = [z(end)-L;z;z(1)+L]; % extend with ghost nodes
    w_ = [w(end);w;w(end)];
    eta = interp1(real(z_),imag(z_),x,'linear',nan);
    phiS = interp1(real(z_),real(w_),x,'linear',nan);
    
    fftEta = fftshift(fft(eta));
    if sum(abs(fftEta(1:floor(end/4))))>.01*sum(abs(fftEta))
        warning('Initial conition may not have been found. Verify that solution exists.')
    end
    
    % phiS0 = ka/k0.*g/omega*sin(xk0-phaseAng);
    % eta0 = ka/k0*(cos(xk0-phaseAng));
    % % phi = ka/k0.*g/omega*sin(xk0-phaseAng)*cosh(k*(h+z))/cosh(k*h);
    % u0 = ka/k0.*g/omega*k0*cos(xk0-phaseAng);
    % v0 =  ka/k0.*g/omega*k0*sin(xk0-phaseAng)*tanh(k0*h);
    % figure('color','w')
    % subplot(311), plot(x,eta,'-',x,eta0,'--');ylabel('\eta'); grid on
    % subplot(312), plot(x,phiS,'-',x,phiS0,'--');ylabel('\phi^S'); grid on
    % subplot(313), plot(x,u0,'-r',x,v0,'-b',real(z),real(dwdz)+out.c_e,'--r',real(z),-imag(dwdz),'--b');ylabel('\phi^S'); grid on
end


if strcmp(surfaceMethod,'decayingConformal')
    
    [eta_adj,H] = initializeInitCond(x,eta,h,10);
    k_cut = k_cut_conformal;
    W = fConformal(x,eta_adj,H,inf);
    phiS_adj = interp1([x-L;x;x+L],[phiS;phiS;phiS],real(W),'linear',nan);
    
%     figure('color','w');
%     f0 = fConformal(x,eta,H,inf);
%     fH = fConformal(x-1i*H,eta_adj,H,inf);
%     -imag(fH(1,:))
%     subplot(2,1,1); plot(x,eta,'-',real(f),imag(f),'--',real(f0),imag(f0),':','linewidth',1.5);ylabel('\eta');title('IC verification')
%     subplot(2,1,2); plot(x,phiS,'-',real(f),phiS_adj,'--','linewidth',1.5);ylabel('\phi^S');
else
    [phiS_adj,eta_adj] = deal(phiS,eta);
    [tInit,yInit] = deal([]);
    H = h;
    k_cut = k_cutTaylor;
end
kx = getKx(x);
% phiS_ip = phiS;
% eta_ip = eta;


tic
[t,y] = ode45(@HOSODE45 ,[t0,t_end],[phiS_adj;eta_adj],ODEoptions);
fprintf('CPU time: %gs\n',toc);

% if t(end) < t_end, return; end

iNaN = find(isnan(y(:,1)),1,'first');
if ~isempty(iNaN), t(iNaN:end)=[]; y(iNaN:end,:)=[]; end
phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);
% interpolate to perscribed times


% t_ip = (0:dt:t_end)';
t_ip = linspace(0,t(end),10).';
nPannel = length(t_ip);
phiS_ip = interp1(t,phiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

if strcmp(surfaceMethod,'decayingConformal')
    W = fConformal(x,eta_ip,H,k_cut);
        
%     fH = fConformal(x-1i*H,eta_ip,H,k_cut);
%     -imag(fH(1,:))
%     figure, plot(x/lambda,(real(fH)-x)/lambda)
    
    for it = i_detailedPlot
        hfCont = figure('color','w'); hold on
        y0 = max(-x(end)/NWaves/2,-H);
        if y0>-H,nan_ = nan; else nan_=1; end
        
        [xi2,sig2]= ndgrid(x,linspace(y0,0,100));
        f2 = fConformal(xi2+1i*sig2,eta_ip(:,it),H,k_cut);
        haCont(1) = subplot(2,1,1);
        contour(real(f2),imag(f2),xi2,30,'r'); hold on
        contour(real(f2),imag(f2),sig2,'b');
        plot([W(:,it);nan;nan_*f2([1,end],1)],'k','linewidth',2)

        haCont(2) = subplot(2,1,2);
        omega = ifft(  fft(phiS_ip(:,it)).*exp(-sig2.*kx).*2./(exp(2*kx.*H)+1).*(abs(kx)<k_cut));
        contourf(real(f2),imag(f2),real(omega),20); hold on
        contour(real(f2),imag(f2),imag(omega),20,'k')
        plot([W(:,it);nan;nan_*f2([1,end],1)],'k','linewidth',2)
        axis(haCont,'equal','off')
        
        if isfinite(H), title(haCont(1), sprintf('H[\\zeta] = %.4g, H[z] = %.4g',H,-imag(f2(1,1)))); end
        
    end

    x_ip = real(W);
%     eta_ip = imag(f); %per definition
else
    x_ip = x+0*eta_ip;
end


[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('%s Tramp%g ka=%.3g,M=%d,H=%.1f',surfaceMethod,Tramp,ka,M,H)),[],[0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),x_ip(:,i),eta_ip(:,i),'k');
    ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),nonLinRamp(t_ip(i))))
    grid(ha(i),'on');
%     axis(ha(i),'equal')
end
% linkaxes(ha)
% ylim(max(res.eta(:))*[-1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)
xlim(ha,x([1,end]));%set(ha,'XLim',x([1,end]));

    
fileName = sprintf('%s%ska%.2g_M%d_h%.2f_Nw%d_dt%.3gT_nx%d_pad%d_ikCut%d',exportPrefix,surfaceMethod,ka,M,h,NWaves,NT_dt,nx,DO_PADDING,round(k_cut/k0)); fileName(fileName=='.')='p';
if DO_EXPORT
    copyfile('./proto_SSGWInit.m',[exportPath,'/script_',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf','-png');
end
if EXPORT_MAT, save([exportPath,'/',fileName]); end




% hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('%s Tramp%g ka=%.3g,M=%d,CPU=%.3g',surfaceMethod,Tramp,ka,M,CPUTime));%[-1587 511 560 1000]
% for iP = 1:nPannel
%     ha(iP) = subplot(nPannel,1,nPannel-iP+1); plot(x_ip(:,iP),eta_ip(:,iP),'k');hold on
%     ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))))
%     box off; grid on;
% end
% linkaxes(ha,'xy');
% set(ha(2:end),'XTick',[]);
% xlim(x([1,end]));