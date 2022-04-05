clear 
clear global
global surfaceMethod timeReached dW DO_PADDING chalikov taylor t_end dim
timeReached = 0; 
addpath c:/gits/timsas2/matlabLibs

g=9.81;
% g=1.0;



%% input
% surfaceMethod = 'Chalikov'; % Chalikov method
surfaceMethod = 'Taylor';  % normal HOS

DO_PADDING = 0;
nx__wave = 2^5;

chalikov.M = 2^11;
RK4dt = 0;%2.5e-3; % set to zero to use ODE45
chalikov.r = .25; % dim.less, =.25 in Chalikov
chalikov.kd__kmax = .5; % .5 in Chalikov


% statistics
rng(500); % gives the seed to all random number generators
waveRandomnessType = 'Rayleigh'; %'Rayleigh' or 'deterministic'



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


h = inf;% 5*lambda; % water depth. 

% waves
% NWaves = 100;
L = 50;
Tp = 1;
Hs = .075;
gamma = 5;
omega0 = 2*pi/Tp;
if isfinite(h)
    k0 = findWaveNumbers(omega0,h,U_curr,0);
else
    if U_curr==0, k0=omega0^2/g; else, k0=(g+2*U_curr*omega0-sqrt(g^2+4*g*U_curr*omega0))/(2*U_curr^2);end
end




% Simulation/plotting time
dt = 200*Tp;
t_end = 9*dt;

    
% stability
INIT_WAVE_TYPE = 'linear';
% packet = 1;
stillWindowWidth__L = 0;%-inf;% 1/4;
stillWindowCentre__L = 2/3;
taperWidth__L = .1*1/4;



taylor.M = 5; 
TRamp = 10*Tp;
taylor.nonLinRamp = @(t) max(0,1-exp(-(t/TRamp)^2));
taylor.k_cut = (taylor.M+5)*k0;
taylor.h = h;

%% Simulation
% L = NWaves*2*pi/k0;
% c_p = 2*pi/T/k0;
NWaves_temp = L*k0/(2*pi);
nx = round(nx__wave*NWaves_temp);
NWaves = nx/nx__wave;
nx = nx__wave*NWaves;
nx = 2*floor(nx/2);


% if strcmp(surfaceMethod,'Chalikov')
%     nx = 2*chalikov.M+1;
% end

dx = L/nx;
x = (0:nx-1)'*dx;
if strcmp(surfaceMethod,'Taylor')
    fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-taylor.k_cut/ ( (2*pi/L)*nx/2),0) )
end
kx = getKx(x);
% stillWindowFilter = 1-exp(-((x-stillWindowCentre__L*L)/stillWindowWidth).^2);
xc = stillWindowCentre__L*L;
tw = stillWindowWidth__L*L/2;
sig = taperWidth__L*L;
stillWindowFilter = 1-.5*(tanh((x-(xc-tw))/sig)-tanh((x-(xc+tw))/sig));

t0 = 0;
initialStepODE = 1e-3*Tp;
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
        nk = floor((nx-1)/2) +1;
        k_plus = (2*pi/L)*(0:nk-1)';
        if isfinite(h)
            f = (k_plus*U_curr+sqrt(g*k_plus.*tanh(k_plus*h)))/(2*pi); 
        else
            f = (k_plus*U_curr+sqrt(g*k_plus))/(2*pi); 
        end
        S = level1.wave.computeWaveSpectralDensity_JONSWAP(f,Hs,Tp,gamma);
        
        df = [f(2)-f(1);.5*(f(3:end)-f(1:end-2)); f(end)-f(end-1)];
        switch waveRandomnessType
            case 'Rayleigh'
                R = (randn(nk,1)+1i*randn(nk,1))/sqrt(2);
            case 'deterministic'
                R = exp(2i*pi*rand(nk,1));
        end
        heta_plus = sqrt(S.*df).* R;  
        heta_plus(k_plus>10*k0)=0; % low-pass filter
%         heta = .5*[heta_plus;zeros(mod(nk-1,2));conj(heta_plus(end:-1:2))];
        heta = .5*[heta_plus;zeros(nx-2*nk+1,1);conj(heta_plus(end:-1:2))];
%         eta = sum( heta.'.*exp( 1i*kx'.*x ),2); % .*exp( -1i*omega*t0)
        eta = ifft(heta)*nx;
        
%         figure,plot(x,eta)
%         F = abs(heta_plus).^2/dx;
%         figure, plot(k_plus,F)
%         S_realization = abs(heta_plus).^2./df;
%         figure, plot(f,S,f,S_realization,'.')
        
%         eta = ka/k0*(cos(xk0-phaseAng));
%         phiS = ka/k0.*g/omega*sin(xk0-phaseAng);

        hphiS = -heta.*1i.*sign(kx).*sqrt(g./abs(kx)); hphiS(1)=0;
        phiS = ifft(hphiS)*nx;
        
    otherwise % input data file assumed
        filePath = ['./IC/',INIT_WAVE_TYPE,'.mat'];
        assert(isfile(filePath));
        load(filePath)
        assert(x0(end)<=x(end))
        eta = [eta0;zeros(nx-length(x0),1)];
        phiS = [phiS0;zeros(nx-length(x0),1)];
%         figure, plot(x,eta,x,phiS,'--')
end
eta = eta.*stillWindowFilter;
phiS = phiS.*stillWindowFilter;


dim.L  = L/(2*pi);
dim.t = sqrt(dim.L/g);
dim.phi = sqrt(dim.L^3*g);
dim.U = dim.L/dim.t;

%% Run simulation
switch surfaceMethod
    case 'Chalikov'
        [eta_adj,chalikov.H] = initializeInitCond(x,eta,h,100);
        FFTeta_adj = fft(eta_adj);
        
        % filtering the initial conditions proved useful when simulating very highly resolved wave in the process of breaking
        filterIC = false;
        if filterIC
            cutFactor = 1e-8;
            FFTeta0 = FFTeta_adj(1); FFTeta_adj(1)=[];
            [~,imax] = max(abs(FFTeta_adj));
            ic = find(abs(FFTeta_adj(imax+1:end))/max(abs(FFTeta_adj)) < cutFactor,1,'first')+imax;
            assert(mod(nx,2)==1); % odd
            FFTeta_adj(ic:nx-ic) = 0; % low-pass filter;
            FFTeta_adj = [FFTeta0;FFTeta_adj];
        end
        
        %     kx = getKx(x);
        %     FFTphiS_adj = -FFTeta_adj.*1i.*sign(kx).*sqrt(g./abs(kx)); FFTphiS_adj(1)=0;
        
        f = fConformal(x,eta_adj,chalikov.H);
        phiS_adj = interp1([x-L;x;x+L],[phiS;phiS;phiS],real(f),'linear',nan);
        FFTphiS_adj = fft(phiS_adj);
        if filterIC
            FFTphiS0 = FFTphiS_adj(1); FFTphiS_adj(1)=[];
            ic = find(abs(FFTphiS_adj)/max(abs(FFTphiS_adj)) < cutFactor,1,'first');
            FFTphiS_adj(ic:nx-ic) = 0; % low-pass filter;
            FFTphiS_adj = [FFTphiS0;FFTphiS_adj];
        end
        
        %     figure('color','w');
        %     f0 = fConformal(x,eta,H,inf);
        %     f = fConformal(x,eta_adj,H,inf);
        %     fH = fConformal(x-1i*H,eta_adj,H,inf);
        %     -imag(fH(1,:))
        %     subplot(2,1,1); plot(x,eta,'-',real(f),imag(f),'--',real(f0),imag(f0),':','linewidth',1.5);ylabel('\eta');title('IC verification')
        %     legend('target in z','actual in z','without adjustment')
        %     subplot(2,1,2); plot(x,phiS,'-',real(f),phiS_adj,'--',real(f0),phiS,':','linewidth',1.5);ylabel('\phi^S');
        %     legend('target in z','actual in z','without adjustment')
        
        %     chalikov.dim = dim;
        ODEoptions.Vectorized = true;
        tic
        switch chalikov.solverSpace
            case 'Fourier'
                if RK4dt~=0
                    [t,y] = RK4(@HOS_Chalikov,[t0,RK4dt,t_end]/dim.t,[FFTphiS_adj/dim.phi;FFTeta_adj/dim.L]);
                else
                    [t,y] = ode45(@HOS_Chalikov ,[t0,t_end]/dim.t,[FFTphiS_adj/dim.phi;FFTeta_adj/dim.L],ODEoptions);
                end
                phiS = ifft(y(:,1:nx),[],2)*dim.phi; eta = ifft(y(:,nx+1:2*nx),[],2)*dim.L;
            case 'physical'
                if RK4dt~=0
                    [t,y] = RK4(@HOS_Chalikov,[t0,RK4dt,t_end]/dim.t,[ifft(FFTphiS_adj)/dim.phi;ifft(FFTeta_adj)/dim.L]);
                else
                    [t,y] = ode45(@HOS_Chalikov ,[t0,t_end]/dim.t,[ifft(FFTphiS_adj)/dim.phi;ifft(FFTeta_adj)/dim.L],ODEoptions);
                end
                phiS = y(:,1:nx)*dim.phi; eta = y(:,nx+1:2*nx)*dim.L;
        end
        fprintf('CPU time: %gs\n',toc);
        t = t*dim.t;
    case 'Taylor'
        tic
        if RK4dt~=0
            [t,y] = RK4(@HOS_Taylor ,[t0,RK4dt,t_end]/dim.t,[phiS/dim.phi;eta/dim.L]);
        else
            [t,y] = ode45(@HOS_Taylor ,[t0,t_end]/dim.t,[phiS/dim.phi;eta/dim.L],ODEoptions);
        end
        fprintf('CPU time: %gs\n',toc);
        phiS = y(:,1:nx)*dim.phi; eta = y(:,nx+1:2*nx)*dim.L;
end
iNaN = find(isnan(phiS(:,1)),1,'first');
if ~isempty(iNaN), t(iNaN:end)=[]; phiS(iNaN:end,:)=[]; eta(iNaN:end,:)=[]; end
clear y

% t_ip = (0:dt:t_end)';
t_ip = linspace(0,t(end),10).';
% t_ip = linspace(.5*t(end),t(end),10).';

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
        
        [xi2,sig2]= ndgrid(x,linspace(y0,0,200));
        f2 = fConformal(xi2+1i*sig2,eta_ip(:,it),H);
        haCont(1) = subplot(2,1,1);
        nLinesX = 50;
        contour(real(f2),imag(f2),xi2,nLinesX,'r'); hold on
        contour(real(f2),imag(f2),sig2,round(-nLinesX*y0/L),'b');
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
        axis(haCont,'tight','equal','off')
        if isfinite(H), title(haCont(1), sprintf('H[\\zeta] = %.4g, H[z] = %.4g',H,-imag(f2(1,1)))); end
    end
else
    x_ip = x+0*eta_ip;
end



[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('%s TRamp%g Hs=%.3g, Tp=%.2g',surfaceMethod,TRamp,Hs,Tp)),[],[0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),x_ip(:,i),eta_ip(:,i),'k');
    if strcmp(surfaceMethod,'Taylor')
        ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),taylor.nonLinRamp(t_ip(i))))
    else
        ylabel(ha(i),sprintf('t = %.2fs',t_ip(i)));
    end
    grid(ha(i),'on');
end
% axis(ha,'equal','tight')
set(ha,'XLim',[min(x_ip(:)),max(x_ip(:))],'YLim',[min(eta_ip(:)),max(eta_ip(:))])
% set(ha,'DataAspectRatio',[1,1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)



iw = round([.75;1]*size(eta,1));
heta = fft(eta(iw(1):iw(2),:).',[],1);
S_mean = mean(abs(2*heta(1:nk,:)).^2./df,2)/df;

hf2 = figure('color','w','position',[1640 164 600 500],'name',sprintf('%s TRamp%g Hs=%.3g, Tp=%.2g',surfaceMethod,TRamp,Hs,Tp));
ha2(1) = subplot(121);ha2(2) = subplot(122);
xlabel(ha2(1),'k');xlabel(ha2(2),'f');
ylabel(ha2(1),'F');ylabel(ha2(2),'S');
hold(ha2(1),'on');hold(ha2(2),'on');
grid(ha2(1),'on');grid(ha2(2),'on');
legend(ha2(1));legend(ha2(2));
for i=[1,nPannel]
    heta = fft(eta_ip(:,i))/nx;
    heta_plus = 2*heta(1:nk);
%     eta = ifft(heta)*nx;
    F = abs(heta_plus).^2/dx;
    S = abs(heta_plus).^2./df;
    assert(max(abs(diff(x_ip(:,i),2)))<1e-6, 'need to interpolate down onto a regular grid');
    hp2(i,1) = plot(ha2(1),k_plus,F,'DisplayName',"t = "+t_ip(i));
    hp2(i,2) = plot(ha2(2),f,F,'DisplayName',"t = "+t_ip(i));
end

if strcmp(surfaceMethod,'Taylor')
    fileName = sprintf('%s%s_Hs%.3g_Tp%.2g_M%d_h%.2f_Nw%d_dt%.3g_nx%d_pad%d_kCut%.4g',exportPrefix,surfaceMethod,Hs,Tp,taylor.M,h,NWaves,dt,nx,DO_PADDING,taylor.k_cut); fileName(fileName=='.')='p';
else
    fileName = sprintf('%s%s_Hs%.3g_Tp%.2g_M%.4g_h%.2f_Nw%d_dt%.3g_nx%d_pad%d',exportPrefix,surfaceMethod,Hs,Tp,chalikov.M,h,NWaves,dt,nx,DO_PADDING); fileName(fileName=='.')='p';
end
if DO_EXPORT
    copyfile('./proto_spectrum.m',[exportPath,'/',fileName,'.m']) 
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
