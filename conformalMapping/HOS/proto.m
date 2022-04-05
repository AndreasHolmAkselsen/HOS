clear 
clear global
global timeReached DO_PADDING taylor t_end map H1 H2 width_xx xi_cut
timeReached = 0; 

g = 9.81;
% g = 1.0;

%% input

xi_cut = 10;

% mapping input
% map.domainType = 'simple';
map.domainType = 'logstip';
% map.domainType = 'double'; 
H2 = .99;
H1 = 1;

width_x = 2.0*max(H1,H2);

% xL = -2.5*max(H1,H2);
% xR = 2.5*max(H1,H2);


% for wave breaking example
DO_PADDING = 1;
RK4dt = 0;%2.5e-3; % set to zero to use ODE45
ka = .05; % linear wave steepness
%  and  NT_dt =  5 /9/T; lambda = 2*pi; g=1;


relTolODE = 1e-4;% 1e-8;
N_SSGW = 2^12; % number of modes in SSGW solution

% Plot & export options
DO_EXPORT = 1;
EXPORT_MAT = 0;
PLOT_MAP = 1;
exportPrefix = 'testConf_';
exportPath = './figures/';
i_detailedPlot = []; %plot contour plots of frame i. Leave empty to skip


% Wave init specification

kH1 = .5;
NWaves = 30;
k0 = kH1/H1;
lambda = 2*pi/k0;
L = lambda*NWaves;
xL = -L/2;
xR = L/2;


% omega = (1+.5*ka^2)*sqrt(g*k0*tanh(k0*H));
% omega = k0*U_curr+sqrt(g*k0*tanh(k0*H)); T = 2*pi/omega;
% omega = sqrt(g*k0*tanh(k0*.5*(H1+H2))); T = 2*pi/omega;
omega = sqrt(g*k0*tanh(k0*(H1))); T = 2*pi/omega;
% T = 1; omega=2*pi/T;
%  if U_curr==0, k0=omega^2/g;else, k0=(g+2*U_curr*omega-sqrt(g^2+4*g*U_curr*omega))/(2*U_curr^2);end;lambda=2*pi/k0;


% Simulation/plotting time
NT_dt =  1;
dt = NT_dt*T;
t_end = 9*dt;
    
% Initial conditions
INIT_WAVE_TYPE = 'linear';  % 'SSGW' or 'linear'
packageWidth__L = .1;  % set to inf if not simulating wave packets
packageCentre__L = -.25;

nx__wave = 2^6;
taylor.M = 3; 
TRamp = 0*T;
taylor.nonLinRamp = @(t) max(0,1-exp(-(t/TRamp)^2));
taylor.k_cut = (taylor.M+5)*k0;


%% Simulation
nx = nx__wave*NWaves;

dx = L/nx;
% x = (0:nx-1)'*dx;
x = linspace(xL,xR,nx+1)';x(end)=[];
fprintf('Fraction of filtered wavespace: %.3g.\n',  max(1-taylor.k_cut/ ( (2*pi/L)*nx/2),0) )
packet = exp(-((x/L-packageCentre__L)/packageWidth__L).^2);

t0 = 0;
initialStepODE = 1e-3*T;
xk0 = k0.*x;
phaseAng = 0*pi/180;
ODEoptions = odeset('RelTol',relTolODE,'InitialStep',initialStepODE);%,'MaxStep',1e-4);



%% Initial conditions:
switch INIT_WAVE_TYPE
    case 'linear'
        h0 = ka/k0*(cos(xk0-phaseAng));
        phiS0 = ka/k0.*g/omega*sin(xk0-phaseAng);
        
    case 'SSGW'
        [z,dwdz,PP] = SSGW(k0*H,ka,N_SSGW);
        
        if isinf(PP(1)), L_scale = 1/k0; else, L_scale = H; end
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
        h0 = interp1(real(z_),imag(z_),x,'linear',nan);
        phiS0 = interp1(real(z_),real(w_),x,'linear',nan);
        fft_h = fftshift(fft(h0));
        if sum(abs(fft_h(1:floor(end/4))))>.01*sum(abs(fft_h))
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
        h0 = [h0;zeros(nx-length(x0),1)];
        varphiS0 = [varphiS0;zeros(nx-length(x0),1)];
%         figure, plot(x,eta,x,phiS,'--')
end
% Adjust in case of wave packets.
h0 = h0.*packet;
phiS0 = phiS0.*packet;

% dim.L  = L/(2*pi);
% dim.t = sqrt(dim.L/g);
% dim.phi = sqrt(dim.L^3*g);
% dim.U = dim.L/dim.t;


%%  map preparation



% Create inverse function through scattered interpolation:
nArrX = 300;
nArrY = 100;
switch map.domainType
    case 'simple'
        map.fz = @fzSimple;
        map.dfz = @dfzSimple;
        yyUpper = (max(hj)+h1)*pi/2;   % from xi->+inf limit
        yyLower = (min(hj)+h2)*pi/2; % from xi->-inf limit
        yyPlotUpper = yyUpper;
        yyPlotLower = 0;     
    case {'logstip','double'}
        if strcmp(map.domainType,'logstip')
            map.fz = @fzStrip;
            map.dfz = @dfzStrip;
        else
            map.fz = @fzDouble;
            map.dfz = @dfzDouble;
            width_xx = fzero( @(Wxx)-2*fzStrip( -Wxx/2 )-width_x,0);
        end
                
        xxL0 = fzero(@(xx) real(map.fz(xx))-width_x,-1); % numerical inverse, upper left corner of z-domain
        assert(isfinite(xxL0))
        dzdzz_min = min(abs(dfzStrip0(linspace(xxL0,-xxL0,21),H1,H2)));
        yyUpper = 1.4*max(h0)./dzdzz_min;
        yyLower = 1.4*min(h0)./dzdzz_min;
%         yyUpper = 1;
%         yyLower = -1;        
        yyPlotUpper = 1;
        yyPlotLower = -pi;   
        map.H = pi;
    otherwise
        error('Domain type''%s'' not recognised.',map.domainType);
end

% dxxExtra = .1*H1;
dxxExtra = .1;
xxL = fzero(@(xx) real(map.fz(xx+1i*yyUpper))-xL,xL )-dxxExtra; % numerical inverse, upper left corner of z-domain
xxR = fzero(@(xx) real(map.fz(xx+1i*yyLower))-xR,xR )+dxxExtra; % numerical inverse, lower left corner of z-domain

% test that the asymptotic treatment for large xi works properly
% figure; zz = (-20:.05:20)'-.5i*pi;
% subplot(121);plot(real(zz),real(map.fz(zz))); hold on; xi_cut = inf; plot(real(zz),real(map.fz(zz)),'--');xi_cut = 10;
% subplot(122);plot(real(zz),imag(map.fz(zz))); hold on; xi_cut = inf; plot(real(zz),imag(map.fz(zz)),'--');xi_cut = 10;

% zz = initializeInitCondNewton(map.fz,@(zz)dfzStrip0(zz,H1,H2),x+1i*h,linspace(xxL,xxR,nx)',20);
% [eta,xi] = initializeInitCond(map.fz,x,h,linspace(xxL,xxR,nx)',x*0,20);
zzArr = linspace(xxL,xxR,nArrX) + 1i*linspace(yyLower,yyUpper,nArrY)';
zArr = map.fz(zzArr);
finvIp = scatteredInterpolant(  real(zArr(:)), imag(zArr(:)) , zzArr(:),'linear','none');
zzS0 = finvIp(x, h0 );

% interpolate onto regulart xi grid.
map.xi = linspace(real(zzS0(1)),real(zzS0(nx)),nx)';
eta0 = interp1(real(zzS0),imag(zzS0),map.xi);
varphiS0 = interp1( [x-L;x;x+L],[phiS0;phiS0;phiS0],real(map.fz(map.xi+1i*eta0)));
% figure, plot(map.fz(zzS)); hold on, plot(x,h,'--'); plot(map.fz(xi+1i*eta),'.')
% figure, plot(x,phiS,real(map.fz(xi+1i*eta)),varphiS,'--');


if PLOT_MAP
    hf_map = figure('color','w','position',[436 63 637 600]); 
    zz = linspace(xxL,xxR,100) + 1i*linspace(yyPlotLower,yyPlotUpper,100)';
    z = map.fz(zz); 
%     etaIp = interp1(real(zzS0),imag(zzS0),real(zz(1,:)),'linear','extrap');
    
    % plot the z-plane
    haz = subplot(211); hold on;
    title('z-plane');xlabel('x');ylabel('i y');box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
%     [~,hcz] = contourf(real(z),imag(z),real(ww),phiLevels,'LineStyle','none');
    
    zPhi = map.fz(linspace(xxL,xxR,10) + 1i*linspace(yyPlotLower,yyPlotUpper,200)');
    zPsi = map.fz(linspace(xxL,xxR,200) + 1i*linspace(yyPlotLower,yyPlotUpper,10)');
    plot(zPhi,'r','linewidth',1); hold on; plot(zPsi.' ,'b')
%     axis equal
    switch map.domainType
        case {'simple','logstip'}
            patch([real(zPsi(1))*[1,1],0,0,real(zPsi(1,end))*[1,1]],[-1.2*max(H1,H2),-H1,-H1,-H2,-H2,-1.2*max(H1,H2)],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
        case 'double'
            patch([real(zPsi(1))*[1,1],-width_x/2*[1,1],width_x/2*[1,1],real(zPsi(1,end))*[1,1]],[-1.1*max(H1,H2),-H2,-H2,-H1,-H1,-H2,-H2,-1.1*max(H1,H2)],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');
    end
    xlabel('x');ylabel('i y');
    axis tight equal
    plot(haz,  [ zArr(1,:),nan, zArr(:,end).',nan,zArr(end,:),nan,zArr(:,1).'],'--k','linewidth',2 )
    plot(x,h0,'k','linewidth',1.5);
    
    % plot the zz-plane
    hazz = subplot(212); hold on    
    title('\zeta-plane'); xlabel('\xi');ylabel('i \sigma'); box off
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])

    contour(real(zz),imag(zz),real(z),'r','linewidth',1);
    contour(real(zz),imag(zz),imag(z),'b','linewidth',1);
    
    plot(hazz,[zzArr(1,1),zzArr(1,end),zzArr(end,end),zzArr(end,1),zzArr(1,1)],'--k','linewidth',2 )
    plot(zzS0,'k','linewidth',1.5);
    axis tight equal

end





%% Run simulation
% tic
% if RK4dt~=0
%     [t,y] = RK4(@HOS_Taylor ,[t0,RK4dt,t_end]/dim.t,[varphiS0/dim.phi;eta0/dim.L]);
% else
%     [t,y] = ode45(@HOS_Taylor ,[t0,t_end]/dim.t,[varphiS0/dim.phi;eta0/dim.L],ODEoptions);
% end
% fprintf('CPU time: %gs\n',toc);
% t = t*dim.t;
% phiS0 = y(:,1:nx)*dim.phi; h0 = y(:,nx+1:2*nx)*dim.L;
tic
if RK4dt~=0
    [t,y] = RK4(@HOS_Taylor ,[t0,RK4dt,t_end],[varphiS0;eta0]);
else
    [t,y] = ode45(@HOS_Taylor ,[t0,t_end],[varphiS0;eta0],ODEoptions);
end
fprintf('CPU time: %gs\n',toc);
varphiS = y(:,1:nx); eta = y(:,nx+1:2*nx);

iNaN = find(isnan(varphiS(:,1)),1,'first');
if ~isempty(iNaN), t(iNaN:end)=[]; varphiS(iNaN:end,:)=[]; eta(iNaN:end,:)=[]; end
clear y



% t_ip = (0:dt:t_end)';
t_ip = linspace(0,t(end),10).';
% t_ip = linspace(.5*t(end),t(end),10).';

nPannel = length(t_ip);
varphiS_ip = interp1(t,varphiS,t_ip).';
eta_ip  = interp1(t,eta ,t_ip).';

zS_ip = map.fz(map.xi+1i*eta_ip);


[hf, ha] = multi_axes(nPannel,1,figure('color','w','position',[1640 164 1081 814],'name',sprintf('Conformal; Tramp%g ka=%.3g',TRamp,ka)),[],[0,0]);
ha = flipud(ha); set([ha(2:end).XAxis],'Visible','off');% if plotting bottom-to-top
hp = 0*t_ip;
% set([ha(1:end-1).XAxis],'Visible','off');% if plotting top-to-bottom
for i=1:nPannel
    hp(i) = plot(ha(i),zS_ip(:,i),'k');
    ylabel(ha(i),sprintf('t = %.2fs\nw_{nl} = %.2f',t_ip(i),taylor.nonLinRamp(t_ip(i))))
    grid(ha(i),'on');
end
% axis(ha,'equal','tight')
set(ha,'XLim',[min(real(zS_ip(:))),max(real(zS_ip(:)))],'YLim',[min(imag(zS_ip(:))),max(imag(zS_ip(:)))])
% set(ha,'DataAspectRatio',[1,1,1])
xlabel(ha(nPannel),'x [m]','fontsize',11)



fileName = sprintf('%s_ka%.2g_M%d_H%.2f_%.2f_Nw%d_dt%.3gT_nx%d_pad%d_kCut%.4g',exportPrefix,ka,taylor.M,H1,H2,NWaves,NT_dt,nx,DO_PADDING,taylor.k_cut); fileName(fileName=='.')='p';

if DO_EXPORT
    copyfile('./proto.m',[exportPath,'/',fileName,'.m']) 
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




function z = fzDouble(zz)
global width_xx
    assert(~any(real(zz)==0,'all'), 'It is assumed that no xi values equal zero.')
    fzL = -fzStrip( -(zz+width_xx/2));
    fzR = +fzStrip( +(zz-width_xx/2));
    z = (fzL-width_x/2).*(real(zz)<0) + (fzR+width_x/2).*(real(zz)>0);  
end

function dzdzz = dfzDouble(zz)
global width_xx
    assert(~any(real(zz)==0,'all'), 'It is assumed that no xi values equal zero.')
    dfzL = -dfzStrip( -(zz+width_xx/2));
    dfzR = +dfzStrip( +(zz-width_xx/2));
    dzdzz = dfzL.*(real(zz)<0) + dfzR.*(real(zz)>0);  
end

function z = fzSimple(zz)
global H1 H2
    assert(H1>H2,'h1>h2 assumed in the ''simple'' configuration.');
    d = H1-H2;
    z = -1i*H2+  2*d/pi*( sqrt(zz/d).*sqrt(1+zz/d)-log(sqrt(zz/d)+sqrt(1+zz/d)) );
end

function dzdzz = dfzSimple(zz)
global H1 H2
    assert(H1>H2,'h1>h2 assumed in the ''simple'' configuration.');
    d = H1-H2;
    dzdzz = 2/pi*sqrt(zz)./sqrt(d+zz);
end

% function z = fzStrip(zz)
% global H1 H2
%     c = h1/h2;
%     lambda = exp(zz+1i*pi); % surface at imag(zz) = 0, bed at imag(zz) = -pi
%     t = sqrt((lambda-c^2)./(lambda-1));
%     z = -1i*h2 + h1/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
% end

function z = fzStrip(zz)
global H1 H2
if H1<H2
    z=fzStrip0(zz,H1,H2);
else
    z=-fzStrip0(-zz,H2,H1);
end
end

function dzdzz = dfzStrip(zz)
global H1 H2
if H1<H2
    dzdzz=dfzStrip0(zz,H1,H2);
else
    dzdzz=-dfzStrip0(-zz,H2,H1);
end
end

function z = fzStrip0(zz,H1,H2)
global xi_cut
    
    iPlus = real(zz) > xi_cut;
    iMinus = real(zz) < -xi_cut;
    iMid = ~(iPlus|iMinus);
    
    z = zeros(size(zz));
    z(iPlus) = fzStrip00(xi_cut,H1,H2) + H2/pi*(zz(iPlus)-xi_cut);
    z(iMinus) = fzStrip00(-xi_cut,H1,H2) + H1/pi*(zz(iMinus)+xi_cut);
    z(iMid) = fzStrip00(zz(iMid),H1,H2);
end

function z = fzStrip00(zz,H1,H2)
    c = H2/H1;
    lambda = -exp(zz); % surface at imag(zz) = 0, bed at imag(zz) = -pi
    t = sqrt((lambda-c^2)./(lambda-1));
    z = -1i*H1 + H2/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
end


function dzdzz = dfzStrip0(zz,H1,H2)
    c2 = (H2/H1)^2; lambda = -exp(zz); 
    t = sqrt((lambda-c2)./(lambda-1));
    %     dzdzz_ =  -lambda*H2/pi.* (c2-1).^2./( (c2-t.^2).*(t.^2-1).*(lambda-1).^2.*sqrt((lambda-c2)./(lambda-1)) );
    dzdzz = H2./(pi*t);
end

function y = log2(x)
    y = log(x);
    ii = imag(x)<0;
    y(ii) = y(ii) + 2i*pi;
end

