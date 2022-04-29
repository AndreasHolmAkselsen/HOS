clear

DO_EXPORT = false;
DO_SUBPLOT = true;
nEv = 400;
axisType = 'normal';%'equal';
PLOT_PRESSURE_INSTEAD = 0;

% HD Wallingford
% h_d = .700+.765;
% h_s = .700;
% Ts = 1./[.25,.5,1];
% exportPath = './stepDoc/Wallingford/';

% Lader
h_d = 1;
h_s = 0.649;
Ts = [1,2,2.5];
exportPath = './stepDoc/Lader/';

% % other
% DO_DEEP_WATER = true;
% h_s = .700;
% Ts = 1./[.25,.5,1];
% exportPath = './stepDoc/other/';

%% code

A = 1;
g = 9.81;

doDeep = exist('DO_DEEP_WATER','var') && DO_DEEP_WATER;
for iT = 1:length(Ts), T=Ts(iT);
    w0 = 2*pi/T;
    t = (0:3)/12 * T;
    % t = (0:3)/8 * T;
    
    if doDeep
        lambda = 2*pi*g/w0^2;
        h_d = 5*lambda;
        x = linspace(-1,1,100)*lambda/2;
        z = linspace(-1,0,100)'*2.5*h_s;
    else
        x = linspace(-1,1,100)*h_d;
        z = linspace(-1,0,100)'*h_d;
    end
	    
    [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s,w0,nEv+1);
    R_n = shiftdim(R_n,-1);
    T_m = shiftdim(T_m,-1);
    k_nv = shiftdim(k_nv,-1);
    k_msv = shiftdim(k_msv,-1);
    k0 = k_nv(1);
    % k_nv.*g.*tanh(k_nv*h_d)-w0^2 == 0
    % k_msv.*g.*tanh(k_msv*h_s)-w0^2 == 0
    
    if PLOT_PRESSURE_INSTEAD
        pf = -1i*w0;
    else
        pf= 1;
    end
    
    ll = lines; quiverColor = ll(2,:);
    
    if DO_SUBPLOT,    hf = figure('color','w','Position',[1668 133 1092 774]); end
    for i_t = 1:length(t)
        phiI = real(pf*A*cosh(k0.*(h_d+z))./cosh(k0.*h_d).*exp(1i*(k0.*x-w0*t(i_t))));
        phiR = real(sum(pf*A*R_n.*cosh(k_nv.*(h_d+z))./cosh(k_nv.*h_d).*exp(1i*(-k_nv.*x-w0*t(i_t))),3));
        phiT = real(sum(pf*A*T_m.*cosh(k_msv.*(h_s+z))./cosh(k_msv.*h_s).*exp(1i*(k_msv.*x-w0*t(i_t))),3));
        phiT0 = real(pf*A*T_m(1).*cosh(k_msv(1).*(h_s+z))./cosh(k_msv(1).*h_s).*exp(1i*(k_msv(1).*x-w0*t(i_t))));
        
        xq = linspace(x(1),x(end),10);
        zq = linspace(z(1),z(end),10)';
        uI = real(1i*k0.*A.*cosh(k0.*(h_d+zq))./cosh(k0.*h_d).*exp(1i*(k0.*xq-w0*t(i_t))));
        uR = real(sum(-1i*k_nv.*A.*R_n.*cosh(k_nv.*(h_d+zq))./cosh(k_nv.*h_d).*exp(1i*(-k_nv.*xq-w0*t(i_t))),3));
        uT = real(sum(1i*k_msv.*A.*T_m.*cosh(k_msv.*(h_s+zq))./cosh(k_msv.*h_s).*exp(1i*(k_msv.*xq-w0*t(i_t))),3));
        uT0 = real(1i*k_msv(1).*A.*T_m(1).*cosh(k_msv(1).*(h_s+zq))./cosh(k_msv(1).*h_s).*exp(1i*(k_msv(1).*xq-w0*t(i_t))));
        wI = real(A*k0.*sinh(k0.*(h_d+zq))./cosh(k0.*h_d).*exp(1i*(k0.*xq-w0*t(i_t))));
        wR = real(sum(A*R_n.*k_nv.*sinh(k_nv.*(h_d+zq))./cosh(k_nv.*h_d).*exp(1i*(-k_nv.*xq-w0*t(i_t))),3));
        wT = real(sum(A*T_m.*k_msv.*sinh(k_msv.*(h_s+zq))./cosh(k_msv.*h_s).*exp(1i*(k_msv.*xq-w0*t(i_t))),3));
        wT0 = real(A*T_m(1).*k_msv(1).*sinh(k_msv(1).*(h_s+zq))./cosh(k_msv(1).*h_s).*exp(1i*(k_msv(1).*xq-w0*t(i_t))));
        
        X = x+0*z;
        phiT(X<=0)=0;
        phiR(X>=0)=0;
        phiI(X>=0)=0;
        phiT0(X<=0)=0;
        
        
        X = xq+0*zq;
        uT(X<=0)=0;
        uR(X>=0)=0;
        uI(X>=0)=0;
        uT0(X<=0)=0;
        wT(X<=0)=0;
        wR(X>=0)=0;
        wI(X>=0)=0;
        wT0(X<=0)=0;
        
        phi = phiI+phiR+phiT;
        u = uI+uR+uT;
        w = wI+wR+wT;
        u(xq>0&zq<-h_s)=nan;
        w(xq>0&zq<-h_s)=nan;
        
        if DO_SUBPLOT
            ha = subplot(ceil(length(t)/2),2,i_t);
        else
            hf = figure('color','w'); ha = gca;
        end
        contourf(x,z,phi);
        hold on;
        quiver(xq,zq,u,w,'color',quiverColor)
        % plot([0,0,x(end)],[-h_d,-h_s,-h_s],'k','linewidth',2)
        patch([0,0,1,1]*x(end),[-h_d,-h_s,-h_s,-h_d],[1,1,1]*.5)
        axis(axisType,[x([1,end]),z([1,end])'])
        
        scale = .4;
        cornerBoxPos = [ha.Position(1:2)+[0,(1-scale)*ha.Position(4)],scale*ha.Position(3:4)];
        has = axes('Position',cornerBoxPos,'box','on');
        contourf(x,z,phiI+phiT0);
        hold on;
        quiver(xq,zq,uI+uT0,wI+wT0,'color',quiverColor)
        patch([0,0,1,1]*x(end),[-h_d,-h_s,-h_s,-h_d],[1,1,1]*.5)
        axis(axisType,[x([1,end]),z([1,end])'],'off')
        plot(has.XLim([1,1,2,2,1]),has.YLim([1,2,2,1,1])','k','linewidth',1.5)
        
        if ~DO_SUBPLOT && DO_EXPORT
            if ~isfolder(exportPath), mkdir(exportPath); end
            if doDeep
                fileName = sprintf('T%.2g_nEv%d_deep_hs%.2g_toT%.2f',T,nEv,h_s,t(i_t)/T); fileName(fileName=='.') = 'p';
            else
                fileName = sprintf('T%.2g_nEv%d_hd%.2g_hs%.2g_toT%.2f',T,nEv,h_d,h_s,t(i_t)/T); fileName(fileName=='.') = 'p';
            end
            savefig(hf,[exportPath,fileName]);
            export_fig(hf,[exportPath,fileName],'-pdf');
        else
            title(sprintf('t = %.2fT',t(i_t)/T))
        end
    end
    if DO_SUBPLOT && DO_EXPORT
        if ~isfolder(exportPath), mkdir(exportPath); end
        if doDeep
            fileName = sprintf('T%.2g_nEv%d_deep_hs%.2f',T,nEv,h_s); fileName(fileName=='.') = 'p';
        else
            fileName = sprintf('T%.2g_nEv%d_hd%.2g_hs%.2f',T,nEv,h_d,h_s); fileName(fileName=='.') = 'p';
        end
        savefig(hf,[exportPath,fileName]);
        export_fig(hf,[exportPath,fileName],'-pdf');
    end
    
end