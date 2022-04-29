clear 
clear global
global x surfaceMethod timeReached dW DO_PADDING kx chalikov taylor
timeReached = 0; 

% g=9.81;
g=1.0;

%% input
surfaceMethod = 'Chalikov'; % Chalikov method
% surfaceMethod = 'Taylor';  % normal HOS

DO_PADDING = 1;
chalikov.M = 3072;
chalikov.dt = 5e-4; % comment out to use ODE45
chalikov.r = .25; % dim.less, =.25 in Chalikov
chalikov.kd__kmax = 0/chalikov.M;%.1; % .5 in Chalikov
ka = .5; % linear wave steepness



chalikov.r = .1; % dim.less, =.25 in Chalikov
chalikov.kd__kmax = -50/chalikov.M;%.1; % .5 in Chalikov
% kd = chalikov.kd__kmax*M;
% k = abs([0:M,-M:-1]');
% mu = chalikov.r*M*((k-kd)/(chalikov.M-kd)).^2.*(k>kd);
%  plot(k,mu)

% DO_PADDING = 0;
% chalikov.M = 3072;
% chalikov.dt = 5e-4;
% chalikov.r = .25; % dim.less, =.25 in Chalikov
% chalikov.kd__kmax = .1; % .5 in Chalikov
% ka = .4; % linear wave steepness


relTolODE = 1e-4;% 1e-8;
N_SSGW = 2^12; % number of modes in SSGW solution

% Plot & export options
DO_EXPORT = 0;
EXPORT_MAT = 0;
PLOT_CURRENT = false;
exportPrefix = '';
exportPath = './figures/';
i_detailedPlot = []; %plot contour plots of frame i. Leave empty to skip




% current specification
U_curr = 0;
currentMatFile = [];
% currentMatFile = '../currendDatabase/uniform_Lh84p8694x2_nx1024nk512_patchPos50xm1p8x15x3p5_Uj0p55_posZ1.mat';
% exportPrefix = 'vortex_';


% NWaves = 60;
NWaves = 1;
% Wave init specification
lambda = 2*pi; k0 = 2*pi/lambda;



h = inf;% 5*lambda; % water depth. 

omega = k0*U_curr+sqrt(g*k0*tanh(k0*h)); T = 2*pi/omega;
% T = 1; omega=2*pi/T;
%  if U_curr==0, k0=omega^2/g;else, k0=(g+2*U_curr*omega-sqrt(g^2+4*g*U_curr*omega))/(2*U_curr^2);end;lambda=2*pi/k0;


% Simulation/plotting time
NT_dt = 10/9/T; %NWaves/9
dt = NT_dt*T;
t_end = 9*dt;

    
% stability
INIT_WAVE_TYPE = 'linear';  % 'SSGW'; 'linear'
packet = 1;
% packageWidth = 3*lambda; packet = exp(-((x-.25/2*L)/packageWidth).^2);
% packageWidth = 3*lambda; packet = exp(-((x/packageWidth-2.5)).^2);
Tramp = 0;%10*T;

% INIT_WAVE_TYPE = SSGW;
% Tramp = 0;


taylor.nx__wave = 2^9;
taylor.M = 3; 
taylor.nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
taylor.k_cut = (taylor.M+5)*k0;

if strcmp(surfaceMethod,'Taylor')
    fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-k_cut/ ( (2*pi/L)*nx/2),0) )
end

%% Simulation

L = NWaves*lambda;
% omega = (1+.5*ka^2)*sqrt(g*k0*tanh(k0*H));
% c_p = 2*pi/T/k0;

if strcmp(surfaceMethod,'Chalikov')
    nx = 2*chalikov.M+1;
else
    nx = taylor.nx__wave*NWaves;
end
dx = L/nx;
x = (0:nx-1)'*dx;

t0 = 0;
initialStepODE = 1e-3*T;
xk0 = k0.*x;
phaseAng = 0*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);%,'MaxStep',1e-4);

%% Specify background current
if ~isempty(currentMatFile)
    load(currentMatFile,'nwt');
    dW = @(zz) nwt.current.dfInterp(real(zz),imag(zz));
else % manual set
    zeta_j = []; A_j = [];
    nMirror = 3; % number of times the domain is repeated in x.

    % single vortex
%     zeta_j = [.5-.075i  ]*L;% object centre
%     A_j    = [ -.2i  ];% object strength
    
    A_j = shiftdim(A_j,-1); zeta_j = shiftdim(zeta_j,-1);% ID_j = shiftdim(ID_j,-1);
    zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);
    
    % % vortex/source/sink
%     A_j = .5*F_j.*c_p.*abs(imag(zeta_j));
    if isempty(A_j)
        W  = @(zz) U_curr.*zz;
        dW = @(zz) U_curr;
    else
        W  = @(zz) sum(A_j.*log(zz-zeta_j) + conj(A_j.*log(conj(zz)-zeta_j)),3:4) + U_curr.*zz;
        dW = @(zz) sum(A_j./(zz-zeta_j) + conj(A_j./(conj(zz)-zeta_j)),3:4) + U_curr;
    end
    
    chalikov.doCurr = ~ (all(isempty(A_j)&&U_curr==0));
end


if PLOT_CURRENT && ~isempty(F_j)
    z = 0:dx:1; z = [-z(end:-1:2),z];
    hf_c = figure('color','w'); ha = gca;
    % plot velocity intensity |U|
    absU = abs(dW(x+1i*z));
    ULim = 2.5*max(abs(dW(x)));
    hIm = imagesc(x',z',absU',[0,ULim]); 
    hIm.AlphaData = (absU<ULim)';
    % plot streamlines
    hold on, axis equal xy
    if exist('W','var'),contour(x',z',imag(W(x+1i*z))',20,'k');end
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
switch INIT_WAVE_TYPE
    case 'linear'
        eta = ka/k0*(cos(xk0-phaseAng));
        phiS = ka/k0.*g/omega*sin(xk0-phaseAng);
        
    case 'SSGW'
        
        [z,dwdz,PP] = SSGW(k0*h,ka,N_SSGW);
        
        if isinf(PP(1)), L_scale = 1/k0; else, L_scale = h; end
        out.c_e = PP(4)*sqrt(g*L_scale); % phase velocity observed from where the meam velocity at the bed is zero
        out.c_s = PP(5)*sqrt(g*L_scale); % mean flow velocity (phase velocity in frame without mean flow)
        out.k = PP(2)/L_scale;
        z = z*L_scale;
        
        % to move wave to centre (optional)
        z = [ z(N_SSGW+1:end)-lambda/2 ; z(1:N_SSGW)+lambda/2 ];
        dwdz = [ dwdz(N_SSGW+1:end); dwdz(1:N_SSGW) ];
        
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
        
        %     phiS0 = ka/k0.*g/omega*sin(xk0-phaseAng);
        %     eta0 = ka/k0*(cos(xk0-phaseAng));
        %     % phi = ka/k0.*g/omega*sin(xk0-phaseAng)*cosh(k*(h+z))/cosh(k*h);
        %     u0 = ka/k0.*g/omega*k0*cos(xk0-phaseAng);
        %     v0 =  ka/k0.*g/omega*k0*sin(xk0-phaseAng)*tanh(k0*h);
        %     figure('color','w')
        %     subplot(311), plot(x,eta,'-',x,eta0,'--');ylabel('\eta'); grid on
        %     subplot(312), plot(x,phiS,'-',x,phiS0,'--');ylabel('\phi^S'); grid on
        %     subplot(313), plot(x,u0,'-r',x,v0,'-b',real(z),real(dwdz)+out.c_e,'--r',real(z),-imag(dwdz),'--b');ylabel('velocity'); grid on
        
    otherwise % input data file assumed
        filePath = ['./IC/',INIT_WAVE_TYPE,'.mat'];
        assert(isfile(filePath));
        load(filePath)
        assert(x0(end)<=x(end))
        eta = [eta0;zeros(nx-length(x0),1)];
        phiS = [phiS0;zeros(nx-length(x0),1)];
%         figure, plot(x,eta,x,phiS,'--')
end
eta = eta.*packet;
phiS = phiS.*packet;


if strcmp(surfaceMethod,'Chalikov')
    
    [eta_adj,chalikov.H] = initializeInitCond(x,eta,h,100);
    
    cutFactor = 1e-8;
    FFTeta_adj = fft(eta_adj); FFTeta0 = FFTeta_adj(1); FFTeta_adj(1)=[];
    ic = find(abs(FFTeta_adj)/max(abs(FFTeta_adj)) < cutFactor,1,'first');
    assert(mod(nx,2)==1); % odd
    FFTeta_adj(ic:nx-ic) = 0; % low-pass filter;
    FFTeta_adj = [FFTeta0;FFTeta_adj];
    
    
%     kx = getKx(x);
%     FFTphiS_adj = -FFTeta_adj.*1i.*sign(kx).*sqrt(g./abs(kx)); FFTphiS_adj(1)=0;

    
    f = fConformal(x,eta_adj,chalikov.H);
    phiS_adj = interp1([x-L;x;x+L],[phiS;phiS;phiS],real(f),'linear',nan);
    
    FFTphiS_adj = fft(phiS_adj); FFTphiS0 = FFTphiS_adj(1); FFTphiS_adj(1)=[];
    ic = find(abs(FFTphiS_adj)/max(abs(FFTphiS_adj)) < cutFactor,1,'first');
    FFTphiS_adj(ic:nx-ic) = 0; % low-pass filter;
    FFTphiS_adj = [FFTphiS0;FFTphiS_adj];
    
%     chalikov.H;
%     FFTeta_adj = zeros(nx,1);
%     FFTeta_adj([2,end]) = .5*ka/k0*nx;
%     FFTphiS_adj = zeros(nx,1);
%     FFTphiS_adj(2) = -.5i*ka/k0.*g/omega*nx;
%     FFTphiS_adj(end)=-FFTphiS_adj(2);
    
    
%     
%     figure('color','w');
%     f0 = fConformal(x,eta,H,inf);
%     f = fConformal(x,eta_adj,H,inf);
%     fH = fConformal(x-1i*H,eta_adj,H,inf);
%     -imag(fH(1,:))
%     subplot(2,1,1); plot(x,eta,'-',real(f),imag(f),'--',real(f0),imag(f0),':','linewidth',1.5);ylabel('\eta');title('IC verification')
%     legend('target in z','actual in z','without adjustment')
%     subplot(2,1,2); plot(x,phiS,'-',real(f),phiS_adj,'--',real(f0),phiS,':','linewidth',1.5);ylabel('\phi^S');
%     legend('target in z','actual in z','without adjustment')
    

    dim.L  = L/(2*pi);
    dim.t = sqrt(dim.L/g);
    dim.phi = sqrt(dim.L^3*g);
    chalikov.t_end = t_end/dim.t;
    chalikov.dim = dim;
    ODEoptions.Vectorized = true;
    tic
    if isfield(chalikov,'dt')
        [t,y] = RK4(@HOSODEeq_mode,[t0/dim.t,chalikov.dt,chalikov.t_end],[FFTphiS_adj/dim.phi;FFTeta_adj/dim.L]);
    else
        [t,y] = ode45(@HOSODEeq_mode ,[t0,t_end]/dim.t,[FFTphiS_adj/dim.phi;FFTeta_adj/dim.L],ODEoptions);
    end
    fprintf('CPU time: %gs\n',toc);
    iNaN = find(isnan(y(:,1)),1,'first');
    if ~isempty(iNaN), t(iNaN:end)=[]; y(iNaN:end,:)=[]; end
    phiS = ifft(y(:,1:nx),[],2)*dim.phi; eta = ifft(y(:,nx+1:2*nx),[],2)*dim.L;
    t = t*dim.t;
else
    kx = getKx(x);
    tic
    [t,y] = ode45(@HOSODE45 ,[t0,t_end],[phiS;eta],ODEoptions);
    fprintf('CPU time: %gs\n',toc);   
    phiS = y(:,1:nx); eta = y(:,nx+1:2*nx);
end
clear y

% t_ip = (0:dt:t_end)';
t_ip = linspace(0,t(end),10).';
% t_ip = linspace(.4*t(end),.7*t(end),10).';

nPannel = length(t_ip);
phiS_ip = interp1(t,phiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

if strcmp(surfaceMethod,'Chalikov')
    W_ip = fConformal(x,eta_ip,chalikov.H);
    x_ip = real(W_ip);
%     eta_ip = imag(f); %per definition
        
%     fH = fConformal(x-1i*H,eta_ip,H,k_cut);
%     -imag(fH(1,:))
%     figure, plot(x/lambda,(real(fH)-x)/lambda)
    
    for it = i_detailedPlot
        H = chalikov.H;
        hfCont = figure('color','w'); hold on
        y0 = max(-x(end)/NWaves/3,-H);
        if y0>-H,nan_ = nan; else nan_=1; end
        
        [xi2,sig2]= ndgrid(x,linspace(y0,0,100));
        f2 = fConformal(xi2+1i*sig2,eta_ip(:,it),H);
        haCont(1) = subplot(2,1,1);
        contour(real(f2),imag(f2),xi2,30,'r'); hold on
        contour(real(f2),imag(f2),sig2,'b');
        plot([W_ip(:,it);nan;nan_*f2([1,end],1)],'k','linewidth',2)

        haCont(2) = subplot(2,1,2);
        kx = getKx(x);
        if isfinite(H)
            omega = ifft(  fft(phiS_ip(:,it)).*exp(-sig2.*kx).*2./(exp(2*kx.*H)+1));
        else
            omega = ifft(  fft(phiS_ip(:,it)).*exp(-sig2.*kx.*(kx<0)).*2.*(kx<0));
        end
        contourf(real(f2),imag(f2),real(omega),20); hold on
        contour(real(f2),imag(f2),imag(omega),20,'k')
        plot([W_ip(:,it);nan;nan_*f2([1,end],1)],'k','linewidth',2)
        axis(haCont,'equal','off')
        if isfinite(H), title(haCont(1), sprintf('H[\\zeta] = %.4g, H[z] = %.4g',H,-imag(f2(1,1)))); end
    end
else
    x_ip = x+0*eta_ip;
end


[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('%s Tramp%g ka=%.3g',surfaceMethod,Tramp,ka)),[],[0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),x_ip(:,i),eta_ip(:,i),'k');
    if strcmp(surfaceMethod,'Taylor')
        ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),nonLinRamp(t_ip(i))))
        fileName = sprintf('%s%s_ka%.2g_M%d_h%.2f_Nw%d_dt%.3gT_nx%d_pad%d_kCut%.4g',exportPrefix,surfaceMethod,ka,taylor.M,h,NWaves,NT_dt,nx,DO_PADDING,taylor.k_cut); fileName(fileName=='.')='p';
    else
        ylabel(ha(i),sprintf('t = %.2fs',t_ip(i)));
        fileName = sprintf('%s%s_ka%.2g_M%.4g_h%.2f_Nw%d_dt%.3gT_nx%d_pad%d',exportPrefix,surfaceMethod,ka,chalikov.M,h,NWaves,NT_dt,nx,DO_PADDING); fileName(fileName=='.')='p';
    end
    grid(ha(i),'on');
end
axis(ha,'equal')
% set(ha,'YLim',[-.65,.65])
% set(ha,'DataAspectRatio',[1,1,1])
linkaxes(ha,'x')
% ylim(max(res.eta(:))*[-1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)
% xlim(ha,[min(x_ip(:)),max(x_ip(:))]);

if strcmp(surfaceMethod,'Taylor')
    fileName = sprintf('%s%s_ka%.2g_M%d_h%.2f_Nw%d_dt%.3gT_nx%d_pad%d_kCut%.4g',exportPrefix,surfaceMethod,ka,taylor.M,h,NWaves,NT_dt,nx,DO_PADDING,taylor.k_cut); fileName(fileName=='.')='p';
else
    fileName = sprintf('%s%s_ka%.2g_M%.4g_h%.2f_Nw%d_dt%.3gT_nx%d_pad%d',exportPrefix,surfaceMethod,ka,chalikov.M,h,NWaves,NT_dt,nx,DO_PADDING); fileName(fileName=='.')='p';
end
if DO_EXPORT
    copyfile('./proto_SSGWInit.m',[exportPath,'/',fileName,'.m']) 
    savefig(hf,[exportPath,'/',fileName]);
    export_fig(hf,[exportPath,'/',fileName],'-pdf','-png');
end
if EXPORT_MAT == 1
    wh = whos;
    vars = setdiff({wh.name},{'t','y','phiS','eta'});
    save([exportPath,'/',fileName],vars{:}); 
elseif EXPORT_MAT == 2
    save([exportPath,'/',fileName]); 
end
