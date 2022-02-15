clear
global h_s h_d Lxx
close all
PLOT_OTHER_PLANES = 0;

h_s = .6;
h_d = 1;


nLines = 10;
nPoints = 100;%500;

Lxx = 10*h_d;


limZeta = [-10*h_d,10*h_d,-pi,.5];
zeta = linspace(limZeta(1),limZeta(2),nPoints) + 1i*linspace(limZeta(3),limZeta(4),nPoints)';


% zeta =  1i*linspace(limZeta(3),limZeta(4),10)';




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
    c = h_d/h_s;
%     c = h_s/h_d;
    t = sqrt((lambda-c^2)./(lambda-1));
    hfT = figure('color','w');hold on
%     contour(real(t),imag(t),real(zeta),'r-','linewidth',1);
%     contour(real(t),imag(t),imag(zeta),'b-','linewidth',1);
    contour(real(t),imag(t),real(lambda),'r-','linewidth',1);
    contour(real(t),imag(t),imag(lambda),'b-','linewidth',1);
    xlabel('Re t');ylabel('i Im t');
    axis equal
    axis([.95,1.05,.0,.05])
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
    box off
end

% plot z-plane
[fz,Lx] = ffz(zeta);
hfz = figure('color','w');hold on

contour(real(fz),imag(fz),real(zeta),nLines,'r','linewidth',1);
contour(real(fz),imag(fz),imag(zeta),nLines,'b','linewidth',1);

hp = patch([min(real(fz(:))).*[1,1],-Lx/2*[1,1],Lx/2*[1,1],max(real(fz(:))).*[1,1]],[-1.1*h_d,-h_d,-h_d,-h_s,-h_s,-h_d,-h_d,-1.1*h_d],.5*[1,1,1],'FaceAlpha',.5,'lineStyle','none');


% map zetaS onto the real equidistant points of the line xS
xS = linspace(min(real(fz(:))), max(real(fz(:))),nLines);
% zeta = linspace(limZeta(1),limZeta(2),100) + 1i*linspace(limZeta(3),limZeta(4),100)';
si = scatteredInterpolant(  real(fz(:)), imag(fz(:)) , zeta(:));
zetaS = si(xS, 0*xS );
plot(real(ffz(zetaS)),imag(ffz(zetaS)),'ok','linewidth',1)
xlabel('x');ylabel('i y');
axis equal tight 
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
box off

hfzz = figure('color','w');
contour(real(zeta),imag(zeta),real(ffz(zeta)),'r','linewidth',1);hold on
contour(real(zeta),imag(zeta),imag(ffz(zeta)),'b','linewidth',1);
plot(zetaS,'-ok','linewidth',1)
axis equal tight 
box off
xlabel('\xi');ylabel('i \sigma');
set(gca,'XAxisLocation','origin','YAxisLocation','origin');%,'XTick',[],'YTick',[])
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



function [fz,width]=ffz(zz)
global Lxx
%     xx0 = real(zz(1,:));
%     zz_minus = zz(:,xx0<0);
%     zz_plus  = zz(:,xx0>0);
%     fzL = -ffz0( -(zz_minus+Lxx/2    ) );
%     fzR = +ffz0( +(zz_plus-Lxx/2) );
%     dx = fzR(1,2)-fzR(1,1);
%     width = 2*real(fzL(1,end)) + dx;
%     fz = [fzL-width/2,fzR+width/2];
    
    
    assert(~any(real(zz)==0,'all'), 'It is assumed that no xi values equal zero.')
    fzL = -ffz0( -(zz+Lxx/2) );
    fzR = +ffz0( +(zz-Lxx/2) );
    width = -2*ffz0( -Lxx/2 );
    fz = (fzL-width/2).*(real(zz)<0) + (fzR+width/2).*(real(zz)>0);  
    
end

function fz=ffz0(zz)
global h_s h_d 
c = h_d/h_s;
lambda = exp(zz+1i*pi); % surface at imag(zz) = 0, bed at imag(zz) = -pi
t = sqrt((lambda-c^2)./(lambda-1));
fz = -1i*h_s + h_d/pi.*(1/c.*log2((t-c)./(t+c))-log((t-1)./(t+1)));
end
