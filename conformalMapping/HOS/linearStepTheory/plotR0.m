clear
g = 9.81;
nEv = 200;

% % for SMB
% h_d = 6;
% h_s = (1:4)';
% T = .5:.1:5;
% exportName = './R0_Lade';

% % for SMB AHA
% h_d = 6;
% h_s = .649*2.4;
% T = .5:.1:5;
% exportName = './R0_SMB';


% % for SMB including short side...
% h_d = 6;
% h_s = .649*[2.4;4]; 
% % h_s = .649*[2.4;3;4;5]; 
% T = .5:.1:5;
% exportName = './R0_SMB_combo';

% % for the Lader tank:
% h_d = 1;
% % h_s = .649;
% h_s = .65;
% T = .6:.05:2.2;


% for SMB including short side...
h_d = 1;
h_s = [.35;.5;.75]; 
% T = 4.1733;
T = linspace(.5,6,100);


w = 2*pi./T;


for i = size(h_s,1):-1:1
    for j = size(T,2):-1:1
%         k0s(i,j) = findWaveNumbers(w0(i,j),h_s(i,j),0,0);
        
         %% [1] linear waves coefficients and wavenumbers
        [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s(i,1),w(1,j),nEv+1);

        R0(i,j) = R_n(1);
        T0(i,j) = T_m(1);
    end
end

hf = figure('color','w'); hf.Position(4) = 400;
grid on; hold on;

% xlabel('Period [s]'); Xvar = T;
% xlabel('$T/\sqrt{H_L/g}$','interpreter','latex'); Xvar = T/sqrt(h_d/g);
xlabel('(kH)_L'); Xvar = findWaveNumbers(w,h_d,0,0)*h_d;
ylabel('R_0');
hp = plot(Xvar,abs(R0),'LineWidth',1);

if length(h_s)>1, legend("H_R/H_L = "+num2str(h_s,'%.2f')+"m",'location','northwest','AutoUpdate','off','fontsize',11); end
% title("h_d = "+h_d+"m");

% % for the Lader tank:
R_shallow = (1-sqrt(h_s./h_d))./(1+sqrt(h_s./h_d));
for R = R_shallow'
    plot(Xvar([1,end]),[1,1].*R,'k--','LineWidth',1);
    text(max(Xvar),R,'Shallow water reflection','fontsize',11,'VerticalAlignment','bottom','HorizontalAlignment','right')
end
xlim([min(Xvar),max(Xvar)])

return

exportName = './R0_k';
% save(exportName)
savefig(hf,exportName)
export_fig(hf,exportName,'-png','-pdf','-m2')


