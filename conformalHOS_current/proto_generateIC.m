clear
global M x k_cut nonLinRamp  surfaceMethod timeReached t_end H dW DO_PADDING kx
timeReached = 0;

%% input
surfaceMethod = 'Chalikov'; % Chalikov method
% surfaceMethod = 'Taylor';  % normal HOS
DO_PADDING = 0;


% Resolution
nx__wave = 2^6;
M = 3; % solution order for Taylor method
relTolODE = 1e-4;% 1e-8;
N_SSGW = 2^12; % number of modes in SSGW solution

% Plot & export options
DO_EXPORT = 1;
EXPORT_MAT = 0;
EXPORT_IC_MAT = 1;
PLOT_CURRENT = false;
exportPrefix = '';
exportPath = './figures/';
i_detailedPlot = []; %plot contour plots of frame i. Leave empty to skip


% % Wave specification
% NWaves = 1;
% lambda = 10;
% ka = .4; % linear wave steepness
% h = 100;%.2*lambda; % water depth. 

% current specification
U_curr = .17;
currentMatFile = [];%'../currendDatabase/workspace_uniform_Lh200x2_nx1024nk512_patchPos100xm1p8x0x0_Uj0_posZ1.mat';



% Wave specification
ka = .2; % linear wave steepness
NWaves = 30;
% lambda = 5; k0 = 2*pi/lambda;
% omega = k0*U_curr+sqrt(g*k0*tanh(k0*h)); T = 2*pi/omega;

T = 1; omega=2*pi/T;
g=9.81; if U_curr==0, k0=omega^2/g;else, k0=(g+2*U_curr*omega-sqrt(g^2+4*g*U_curr*omega))/(2*U_curr^2);end;lambda=2*pi/k0;


% Domain
h = 100;%.2*lambda; % water depth. 

% some computations...
L = NWaves*lambda;
g = 9.81;
% omega = (1+.5*ka^2)*sqrt(g*k0*tanh(k0*H));
c_p = 2*pi/T/k0;
nx = nx__wave*NWaves;
dx = L/nx;
x = (0:nx-1)'*dx;% - L/2;

% Simulation/plotting time
NT_dt = 30/9;
dt = NT_dt*T;
t_end = 9*dt;

    
% stability
INIT_WAVE_TYPE = 'SSGW';  % 'SSGW'; 'linear'
% packet = 1;
packageWidth = 3*lambda; packet = exp(-((x-.5*L)/packageWidth).^2);
Tramp = 10*T;

% INIT_WAVE_TYPE = SSGW;
% Tramp = 0;

nonLinRamp = @(t) max(0,1-exp(-(t/Tramp)^2));
k_cutTaylor = (M+5)*k0;
% k_cutTaylor = 10*(2*pi/L);
% k_cutTaylor = inf;
k_cut_conformal = .25; % dim.less, =.25 in Chalikov

if strcmp(surfaceMethod,'Taylor')
    fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-k_cutTaylor/ ( (2*pi/L)*nx/2),0) )
end
%% Simulation
t0 = 0;
initialStepODE = 1e-3*T;
xk0 = k0.*x;
phaseAng = 0*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);

%% Specify background current
if ~isempty(currentMatFile)
    load(currentMatFile,'nwt');
    dW = @(zz) nwt.current.dfInterp(real(zz),imag(zz));
else % manual set
    zeta_j = []; F_j = [];
    nMirror = 3; % number of times the domain is repeated in x.

    % single vortex
%     zeta_j = [.5-.075i  ]*L;% object centre
%     F_j    = [ -.2i  ];% object strength
    
    F_j = shiftdim(F_j,-1); zeta_j = shiftdim(zeta_j,-1);% ID_j = shiftdim(ID_j,-1);
    zeta_j = zeta_j + L*shiftdim(-nMirror:nMirror,-2);
    
    % % vortex/source/sink
    A_j = .5*F_j.*c_p.*abs(imag(zeta_j));
    if isempty(F_j)
        W  = @(zz) U_curr.*zz;
        dW = @(zz) U_curr;
    else
        W  = @(zz) sum(A_j.*log(zz-zeta_j) + conj(A_j.*log(conj(zz)-zeta_j)),3:4) + U_curr.*zz;
        dW = @(zz) sum(A_j./(zz-zeta_j) + conj(A_j./(conj(zz)-zeta_j)),3:4) + U_curr;
    end
    % doublet
    % A_j = .5*F_j.*c_p.*abs(imag(zeta_j)).^2;
    % f = @(zz) sum(-A_j./(zz-zeta_j) - conj(A_j)./(zz-conj(zeta_j)),3:4);
    % dW = @(zz) sum(A_j./(zz-zeta_j).^2 + conj(A_j)./(zz-conj(zeta_j)).^2,3:4);
    
    if all(isempty(F_j)&&U_curr==0), dW=@(zz)0; end
end


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
        figure, plot(x,eta,x,phiS,'--')
end
eta = eta.*packet;
phiS = phiS.*packet;

if strcmp(surfaceMethod,'Chalikov')
    
    [eta_adj,H] = initializeInitCond(x,eta,h,10);
    k_cut = k_cut_conformal;
    W = fConformal(x,eta_adj,H,inf);
    phiS_adj = interp1([x-L;x;x+L],[phiS;phiS;phiS],real(W),'linear',nan);
    
%     figure('color','w');
%     f0 = fConformal(x,eta,H,inf);
%     f = fConformal(x,eta_adj,H,inf);
%     fH = fConformal(x-1i*H,eta_adj,H,inf);
%     -imag(fH(1,:))
%     subplot(2,1,1); plot(x,eta,'-',real(f),imag(f),'--',real(f0),imag(f0),':','linewidth',1.5);ylabel('\eta');title('IC verification')
%     legend('target in z','actual in z','without adjustment')
%     subplot(2,1,2); plot(x,phiS,'-',real(f),phiS_adj,'--',real(f0),phiS,':','linewidth',1.5);ylabel('\phi^S');
%     legend('target in z','actual in z','without adjustment')
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
% t_ip = linspace(.9*t(end),t(end),10).';

nPannel = length(t_ip);
phiS_ip = interp1(t,phiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

if strcmp(surfaceMethod,'Chalikov')
    W = fConformal(x,eta_ip,H,k_cut);
    x_ip = real(W);
%     eta_ip = imag(f); %per definition
        
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
% return

if EXPORT_IC_MAT
    % % moving and tapering around wave packet.
    for i = 1:nPannel
        [~,imax] = max(abs(eta_ip(:,i)));
        [eta1,phiS1,x1]=deal(eta_ip(:,i),phiS_ip(:,i),x);
        if imax<.25*nx || imax > .75*nx
            eta1 = [eta1(nx/2:end);eta1;eta1(1:nx/2)];
            phiS1 = [phiS1(nx/2:end);phiS1;phiS1(1:nx/2)];
            x1 =  [x(nx/2:end)-L;x(:);x(1:nx/2)+L];
            [~,imax] = max(abs(eta1).*[zeros(nx/2+1,1);ones(nx,1);zeros(nx/2,1)]);
        end
        imax = imax - round(2.5/dx); % adjust!
        tw = 2*packageWidth;
        taper = .5*(tanh(x1-(x1(imax)-tw))-tanh(x1-(x1(imax)+tw)));
%         plot(x1,eta1,x1,phiS1,'--',x1,taper,':')
        
        eta0 = eta1.*taper;
        phiS0 = phiS1-mean(phiS1);
        phiS0 = phiS0.*taper;
        window = [-1,1]*round(3*packageWidth/dx);
        ii = (window(1):window(2))+imax;
        x0 = x1(ii);
        eta0 = eta0(ii);
        phiS0 = phiS0(ii);
%         plot(x0,eta0,x0,phiS0,'--')
        save(['./IC/',fileName,'_',num2str(i)],'x0','eta0','phiS0');
    end
end
% eta0 = ifft(fft(eta0).*(abs(kx)>4*2*pi/L));% high-pass filter
% phiS0 = ifft(fft(phiS0).*(abs(kx)>4*2*pi/L));




 
% hf = figure('color','w','Position',[527  0  1056  1000],'name',sprintf('%s Tramp%g ka=%.3g,M=%d,CPU=%.3g',surfaceMethod,Tramp,ka,M,CPUTime));%[-1587 511 560 1000]
% for iP = 1:nPannel
%     ha(iP) = subplot(nPannel,1,nPannel-iP+1); plot(x_ip(:,iP),eta_ip(:,iP),'k');hold on
%     ylabel(sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(iP),nonLinRamp(t_ip(iP))))
%     box off; grid on;
% end
% linkaxes(ha,'xy');
% set(ha(2:end),'XTick',[]);
% xlim(x([1,end]));