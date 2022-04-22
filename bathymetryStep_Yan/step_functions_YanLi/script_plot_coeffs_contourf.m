clear all
load('R_10.mat')
load('R_20.mat')
load('T_10.mat')
load('T_20.mat')
load('k_0sm.mat')

%%
addpath('C:\Users\yal\OneDrive - Nexus365\Matlab\Functions')



%%
% g        =  9.81;
% nhd      =  512;
% h_d      =  1;
% k_0h_d   =  linspace(0.3, 2*pi, nhd);
% h_ratio  =  linspace(0.2,1,nhd);

%% dimensional parameters
% omega_0  = sqrt(k_0h_d.*tanh(k_0h_d)*g/h_d);
% h_s      = h_d*h_ratio;
% k_0      = k_0h_d/h_d;

%% SUB-harmonic waves here

[h_sm,k_0m]       = meshgrid(h_s,k_0);
omega_0m          = sqrt(k_0m*h_d.*tanh(k_0m*h_d)*g/h_d);
c_g0m             = 0.5*omega_0m./k_0m.*(1+2*k_0m*h_d./sinh(2*k_0m*h_d));
term1B            = (2*g*h_d-c_g0m.^2 )/2./sinh(2*k_0m*h_d) + 2*g.*c_g0m./omega_0m;
B_d               = -1 ./ (4*(g*h_d-c_g0m.^2)) .*term1B; 

k0hs              = h_sm.*k_0sm;
omega_test        = sqrt(k0hs.*tanh(k0hs)*g./h_sm);
c_g0sm            = 0.5*omega_0m./k_0sm.*(1+2*k0hs./sinh(2*k0hs));
term1Bs           = (2*g*h_sm-c_g0sm.^2 )/2./sinh(2*k0hs) + 2*g.*c_g0sm./omega_0m;
cg_longs         =  sqrt(g*h_s) ;
diff_cg           = cg_longs.^2 - c_g0sm.^2 ;
B_s               = -1 ./ (4*(g*h_s-c_g0sm.^2)) .*term1Bs; 

c_gsubR           = sqrt(g*h_d);
c_gsubT           = sqrt(g*h_sm);
gamma_h           = h_sm/h_d;
F_1               = k_0sm.*B_s.*h_sm.*g.*(abs(T_10)).^2./c_g0sm - (1-(abs(R_10)).^2).*B_d.*h_d.*k_0m*g./c_g0m;  % from the momentum
F_2               = k_0sm.*B_s.*(abs(T_10)).^2 - (1+(abs(R_10)).^2).*B_d.*k_0m; 
B_Rf              = (F_1-sqrt(g*h_sm).*F_2)./(-(sqrt(g*h_d)+sqrt(g*h_sm)).*k_0m);
B_Tf              = (F_1+sqrt(g*h_d).*F_2)./(-(sqrt(g*h_d)+sqrt(g*h_sm)).*k_0sm);

B_Rf(isinf(B_Rf)) = 0;
B_Tf(isinf(B_Tf)) = 0;

%%
B_Rf = B_Rf.';
B_Tf = B_Tf.';
B_d = B_d.';
B_s = B_s.';

%%
R_10 = R_10.';
R_20 = R_20.';
T_10 = T_10.';
T_20 = T_20.';
k_0sm = k_0sm.';

ang_R10  = atan(imag(R_10)./real(R_10)) +...
                heaviside(-real(R_10)).*sign(imag(R_10))*pi;
ang_T10  = atan(imag(T_10)./real(T_10)) +...
                heaviside(-real(T_10)).*sign(imag(T_10))*pi;       
ang_R20  = atan(imag(R_20)./real(R_20)) +...
                heaviside(-real(R_20)).*sign(imag(R_20))*pi;
ang_T20  = atan(imag(T_20)./real(T_20)) +...
                heaviside(-real(T_20)).*sign(imag(T_20))*pi;

T_10(isnan(T_10)) = 1;
R_10(isnan(R_10)) = 0;
R_20(isnan(R_20)) = 0;
T_20(isnan(T_20)) = 0;

%% The linear coefficients
addpath('C:\Users\yal\OneDrive - Nexus365\Matlab\Functions')
ymax = 1;
gap = 0.025; % [vertical,horizontal]
marg_h = [0.14 0.12];     % [lower uppper] 
marg_w = [0.042 0.01];      %  [left right] 
% marg_h = 0.142;
% marg_w = 0.041;
omega_0  = sqrt(k_0h_d.*tanh(k_0h_d)*g/h_d);
f1 = figure(1);

f1.Units = 'centimeters';
nwith = 28;
nheigt = 7.5;
% Axis  'tight'
set(f1,'Position',[10 5 nwith nheigt]);
n_row = 1;
n_colum = 4;

i=1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [0 0.05 0.1 0.15 0.2 0.3];
[C,h] = contourf(k_0h_d,h_ratio,abs(R_10),v); colorbar
clabel(C,h,v,'FontName','Times New roman','fontsize',12)
xlabel('$k_0h_d$','Interpreter','latex')
ylabel('$h_s/h_d$','Interpreter','latex')
title('(a) $|R_{0}|$','Interpreter','latex')
set(gca,'fontsize',12,'FontName','Times New roman')
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'fontsize',12)

i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = -[0 0.1 0.2 0.3 0.4 0.5 0.6];
[C,h] = contourf(k_0h_d,h_ratio,ang_R10/pi,v); colorbar
clabel(C,h,v,'FontName','Times New roman','fontsize',12)
xlabel('$k_0h_d$','Interpreter','latex')
title('(b) arg $(R_{0})/\pi$','Interpreter','latex')
set(gca,'fontsize',10,'FontName','Times New roman')
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'fontsize',12)

%
i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [0.9 0.95 1 1.05 1.1 1.2 1.3];
[C,h] = contourf(k_0h_d,h_ratio,abs(T_10),v);colorbar
clabel(C,h,v,'FontName','Times New roman','fontsize',12)
xlabel('$k_0h_d$','Interpreter','latex')
title('(c) $|T_{0}|$','Interpreter','latex')
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'fontsize',12,'FontName','Times New roman')
% set(gca,'ylim',[0.2 ymax],'fontsize',10)

%
i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [0 0.01 0.02 0.03 0.04 0.05 ];
[C,h] = contourf(k_0h_d,h_ratio,ang_T10/pi,v);colorbar
clabel(C,h,v,'FontName','Times New roman','fontsize',12)
xlabel('$k_0h_d$','Interpreter','latex')
title('(d) arg $(T_{0})/\pi$','Interpreter','latex')
set(gca,'fontsize',12,'FontName','Times New roman')
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'fontsize',12)

clear newColormap
intensity = 60; %choose from 1 to 64, where 64 will range to a full black color and 1 will make everything white
newColormap = colormap(flipud(colormap('bone')));
newColormap = newColormap(15:8:intensity,:);
colormap(newColormap);

set(gcf,'color','w');
% export_fig('nonDcoeff_a-d.fig')
% export_fig('nonDcoeff_a-d.pdf')

%%
ymax = 1;
gap = [0.1,0.025]; % [vertical,horizontal]
marg_h = [0.1 0.05];     % [lower uppper] 
marg_w = [0.05 0.005];      %  [left right] 
% marg_h = 0.142;
% marg_w = 0.041;
omega_0  = sqrt(k_0h_d.*tanh(k_0h_d)*g/h_d);
f1 = figure(2);

f1.Units = 'centimeters';
nwith = 28;
nheigt = 14;
% Axis  'tight'
set(f1,'Position',[10 5 nwith nheigt]);
omega_0  = sqrt(k_0h_d.*tanh(k_0h_d)*g/h_d);
n_row = 2;
n_colum = 4;
%
i = 1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [0 0.1 0.2 0.5 0.8 1 1.5];
[C,h] = contourf(k_0h_d(1:end-1),h_ratio,abs(R_20(:,1:end-1))./omega_0(1:end-1),v);colorbar
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'zlim',[0 1.5],'ylim',[0.2 ymax],'fontsize',12)
val_1 = [0 0.1 0.2 0.5 0.8 ]; 
clabel(C,h,val_1,'FontName','Times New roman','fontsize',12)
ylabel('$h_s/h_d$','Interpreter','latex')
title('(e) $|R_{20}|$','Interpreter','latex')

set(gca,'fontsize',12,'FontName','Times New roman')
% keyboard
%
i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [0 1  3 5 10 20];
val_1 = [0 0.1 0.2 0.5 0.8 ]; 
[C,h] = contourf(k_0h_d(1:end-1),h_ratio,abs(T_20(:,1:end-1))./omega_0(1:end-1),v); colorbar
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'zlim',[0 10],'fontsize',12)
clabel(C,h,v,'FontName','Times New roman','fontsize',12)

set(gca,'fontsize',12,'FontName','Times New roman')
title('(f) $|T_{20}|$','Interpreter','latex')
% 'xlim',[0.2 1]*2*pi,

%
%
i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [-5 -3 -2 -1 -0.5  0 1 2 3 5];
val_1 = [ -2 -1 -0.5  0]; 
[C,h] = contourf(k_0h_d(1:end-1),h_ratio,B_Rf(:,1:end-1),v); colorbar
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'zlim',[-5 5],'fontsize',12)
clabel(C,h,val_1,'FontName','Times New roman','fontsize',12)
% set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'zlim',[-3 0],'fontsize',12)

set(gca,'fontsize',12,'FontName','Times New roman')
title('(g) $B_{R}^f$','Interpreter','latex')

i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = -[5 3 2 1 0.5  0];
[C,h] = contourf(k_0h_d(1:end-1),h_ratio,B_d(:,1:end-1),v); colorbar
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'zlim',[-3 0],'fontsize',12)
clabel(C,h,v,'FontName','Times New roman','fontsize',12)

% set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'fontsize',12)
set(gca,'fontsize',12,'FontName','Times New roman')
title('(h) $B_d$','Interpreter','latex')


 i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [0 0.01 0.02 0.05 0.1];
val_1 = [0.01  0.05 0.1];
ang_R20(isnan(ang_R20)) = 0;
[C,h] = contourf(k_0h_d,h_ratio,ang_R20/pi,v);colorbar

clabel(C,h,val_1,'FontName','Times New roman','fontsize',12)
ylabel('$h_s/h_d$','Interpreter','latex')
title('(i) arg $(R_{20})/\pi$','Interpreter','latex')
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'zlim',[0 1.5],'ylim',[0.2 ymax],'fontsize',12)
set(gca,'fontsize',12,'FontName','Times New roman')
xlabel('$k_0h_d$','Interpreter','latex')

%
i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
 v     = [-1 -0.99 -0.98 -0.97 -0.95 -0.9];
 ang_T20(isnan(ang_T20)) = 0;
val_1 =  [-0.99 -0.98 -0.97 -0.95 -0.9];
 ang_T20(ang_T20==0) = -pi;
[C,h] = contourf(k_0h_d,h_ratio,ang_T20/pi,v);colorbar
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'zlim',[0 1.5],'ylim',[0.2 ymax],'fontsize',12)
clabel(C,h,val_1,'FontName','Times New roman','fontsize',12)
title('(j) arg $(T_{20})/\pi$','Interpreter','latex')
set(gca,'fontsize',12,'FontName','Times New roman')
xlabel('$k_0h_d$','Interpreter','latex')

%
i = i+1;
h(i)  = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [ 0 1  3 5 10];
val_1 = [ 0 1  3 5 10];
[C,h] = contourf(k_0h_d(1:end-1),h_ratio,B_Tf(:,1:end-1),v); colorbar
xlabel('$k_0h_d$','Interpreter','latex')
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'zlim',[-8 5],'fontsize',12)
clabel(C,h,val_1,'FontName','Times New roman','fontsize',12) 

% set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'fontsize',12)
set(gca,'fontsize',12,'FontName','Times New roman')
title('(k) $B_{T}^f$','Interpreter','latex')
% 


%
i = i+1;
h(i) = subtightplot(n_row,n_colum,i,gap,marg_h,marg_w);
v     = [-20 -10 -5 -3 -2 -1 -0.5 0];
val_1 = [-20 -10 -5 -3 -2 -1 -0.5 0];
[C,h] = contourf(k_0h_d(1:end-1),h_ratio,B_s(:,1:end-1),v); colorbar
set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'zlim',[-5 5],'fontsize',12)
clabel(C,h,val_1,'FontName','Times New roman','fontsize',12)
% set(gca,'YScale', 'log','XScale', 'log','xtick',[0.5 1 2 4 6],'ylim',[0.2 ymax],'zlim',[-3 0],'fontsize',12)
xlabel('$k_0h_d$','Interpreter','latex')
set(gca,'fontsize',12,'FontName','Times New roman')
title('(l) $B_{s}$','Interpreter','latex')

clear newColormap
intensity = 60; %choose from 1 to 64, where 64 will range to a full black color and 1 will make everything white
newColormap = colormap(flipud(colormap('bone')));
newColormap = newColormap(15:8:intensity,:);
colormap(newColormap);
%%
set(gcf,'color','w');

export_fig('nonDcoeff_e-j.fig')
export_fig('nonDcoeff_e-j.pdf')

%%
% import matplotlib.pyplot as plt
% left = 0.125  # the left side of the subtightplots of the figure
% right = 0.9   # the right side of the subtightplots of the figure
% bottom = 0.1  # the bottom of the subtightplots of the figure
% top = 0.9     # the top of the subtightplots of the figure
% wspace = 0.2  # the amount of width reserved for space between subtightplots,
%               # expressed as a fraction of the average axis width
% hspace = 0.2  # the amount of height reserved for space between subtightplots,
%               # expressed as a fraction of the average axis height
% left = 0.125 ;
% bottom = 0.1  ;
% top = 0.9     ;
% wspace = 0.2  ;              
% hspace = 0.2  ;
              
%ply.subtightplots_adjust(left=0.1, bottom='None', right= 'None', top='None', wspace='None', hspace='None')
% 
% nwith = 24/3;
% nheigt = nwith*2/3;
% nc = 2;
% nr = 3;
% 
% nwith = 24;
% nheigt = nwith*4/5;
% 
% gap = 0.2;
% 
% % for i = 1:nr
% %   for  j = 1:nc
% %      set( h(j+nc*(i-1)), 'Position', [(j-1)*nwith/nc +gap*(j) (i-1)*nwith/nc +gap*(i) 0 nwith/nc-gap*(nc+1) nheigt/nr-gap*(nr+1) ])
% %   end    
% % end
% 
% % set(gca, 'LooseInset', get(gca,'TightInset'))
% saveas(f1,'nonD_coeffs')
% saveas(f1,'nonD_coeffs','epsc')
% saveas(f1,'nonD_coeffs','tif')



