clear
global h1 h2

PLOT_OTHER_PLANES = true;


h = 1;
d = .4;

h1 = h-d;
h2 = h;
nLines = 10;
nPoints = 100;%500;




limZeta = [-5*h,5*h,-pi,.5];
zeta = linspace(limZeta(1),limZeta(2),nPoints) + 1i*linspace(limZeta(3),limZeta(4),nPoints)';


% zeta =  1i*linspace(limZeta(3),limZeta(4),10)';


% plot zeeee-plane
% lsZeta = logspace( log10(min(abs(lambda(:)))),log10(max(abs(lambda(:)))),5);
% contour(real(fz),imag(fz),real(lambda), [-lsZeta,lsZeta]  ,'k-');
% contour(real(fz),imag(fz),imag(lambda), lsZeta,'k-');

% lsT = logspace( log10(min(abs(t(:)))),log10(max(abs(t(:)))),5);
% contour(real(fz),imag(fz),real(t), lsT  ,'g-');
% contour(real(fz),imag(fz),imag(t), lsT,'g-');

if PLOT_OTHER_PLANES
    % Plot lambda plane
    hfLam = figure('color','w');hold on
    lambda = exp(zeta+1i*pi); 
    contour(real(lambda),imag(lambda),real(zeta),'r-','linewidth',1);
    contour(real(lambda),imag(lambda),imag(zeta),'b-','linewidth',1);
    xlabel('Re \lambda');ylabel('i Im \lambda');
    axis equal tight 
    set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
    box off

% plot t-plane
    c = h2/h1;
    t = sqrt((lambda-c^2)./(lambda-1));
    hfT = figure('color','w');hold on
%     contour(real(t),imag(t),real(zeta),'r-','linewidth',1);
%     contour(real(t),imag(t),imag(zeta),'b-','linewidth',1);
    contour(real(t),imag(t),real(lambda),'r-','linewidth',1);
    contour(real(t),imag(t),imag(lambda),'b-','linewidth',1);
    xlabel('Re t');ylabel('i Im t');
    axis equal
    axis([.95,1.05,.0,.05])
    set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
    box off
end

% plot z-plane
fz = ffz(zeta);
hfz = figure('color','w');hold on

contour(real(fz),imag(fz),real(zeta),nLines,'r','linewidth',1);
contour(real(fz),imag(fz),imag(zeta),nLines,'b','linewidth',1);

hp = patch([min(real(fz(:))).*[1,1],0,0,max(real(fz(:))).*[1,1]],[-1.1*h2,-h1,-h1,-h2,-h2,-1.1*h2],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');


% map zetaS onto the real equidistant points of the line xS
xS = linspace(min(real(fz(:))), max(real(fz(:))),nLines);
% zeta = linspace(limZeta(1),limZeta(2),100) + 1i*linspace(limZeta(3),limZeta(4),100)';
si = scatteredInterpolant(  real(fz(:)), imag(fz(:)) , zeta(:));
zetaS = si(xS, 0*xS );
plot(real(ffz(zetaS)),imag(ffz(zetaS)),'ok','linewidth',1)
xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
box off

hfzz = figure('color','w');
contour(real(zeta),imag(zeta),real(ffz(zeta)),'r','linewidth',1);hold on
contour(real(zeta),imag(zeta),imag(ffz(zeta)),'b','linewidth',1);
plot(zetaS,'-ok','linewidth',1)
axis equal tight 
box off
xlabel('\xi');ylabel('i \sigma');
set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[])
return
export_fig(hfz,'CCMei_step','-pdf');
export_fig(hfzz,'CCMei_step_inv','-pdf');
export_fig(hfLam,'CCMei_step_lambda','-pdf');
export_fig(hfT,'CCMei_step_t','-pdf');

% function y = sqrt2(x)
%     y = sqrt(x);
%     ii = imag(x)<0;
%     y(ii) = -y(ii);
% end

function y = log2(x)
    y = log(x);
    ii = imag(x)<0;
    y(ii) = y(ii) + 2i*pi;
end

% function fz=ffz(zeta)
% global h1 h2
% c = h2/h1;
% % lambda = exp(zeta); % surface at imag(zeta) = pi, bed at imag(zeta) = 0
% lambda = -exp(zeta); % surface at imag(zeta) = 0, bed at imag(zeta) = -pi
% t = sqrt((lambda-c^2)./(lambda-1));
% % t = sqrt2(lambda-c^2)./sqrt2(lambda-1);
% fz = -1i*h1 + h2/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
% end



function fz=ffz(zeta)
global h1 h2
c = h2/h1;
lambda = exp(zeta); % surface at imag(zeta) = 0, bed at imag(zeta) = -pi
t = sqrt((lambda+c^2)./(lambda+1));
fz = -1i*h1 + h2/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
end


% A = rand(2000);
% b = rand(2000,1);

% tic
% c=A\b;
% toc

