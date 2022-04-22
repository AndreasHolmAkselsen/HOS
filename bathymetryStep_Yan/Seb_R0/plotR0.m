clear
g = 9.81;
nEv = 100;

% for SMB
h_d = 6;
h_s = (1:4)';
T = .5:.1:5;
exportName = './R0_Lade';

% for SMB AHA
h_d = 6;
h_s = .649*2.4;
T = .5:.1:5;
exportName = './R0_SMB';


% for SMB including short side...
h_d = 6;
h_s = .649*[2.4;4]; 
% h_s = .649*[2.4;3;4;5]; 
T = .5:.1:5;
exportName = './R0_SMB_combo';


% % for the Lader tank:
% h_d = 1;
% % h_s = .649;
% h_s = .65;
% T = .6:.05:2.2;

w0 = 2*pi./T;


for i = size(h_s,1):-1:1
    for j = size(T,2):-1:1
%         k0s(i,j) = findWaveNumbers(w0(i,j),h_s(i,j),0,0);
        
         %% [1] linear waves coefficients and wavenumbers
        [R_n,T_m,k_nv,k_msv] = monochramonic_coefficient_final(h_d,h_s(i,1),w0(1,j),nEv+1);
        
        R0(i,j) = R_n(1);
        T0(i,j) = T_m(1);
    end
end


hf = figure('color','w'); hf.Position(4) = 400;
hp = plot(T,abs(R0),'LineWidth',1);
grid on;
xlabel('Period [s]'); ylabel('R_0');
if length(h_s)>1, legend("D = "+num2str(h_s,'%.2f')+"m",'location','northwest','AutoUpdate','off','fontsize',11); end
% title("h_d = "+h_d+"m");

% % for the Lader tank:
hold on
R_shallow = (1-sqrt(h_s./h_d))./(1+sqrt(h_s./h_d));
for R = R_shallow'
    plot(T([1,end]),[1,1].*R,'k--','LineWidth',1);
    text(max(T),R,'Shallow water reflection','fontsize',11,'VerticalAlignment','bottom','HorizontalAlignment','right')
end

return
save(exportName)
savefig(hf,exportName)
export_fig(hf,exportName,'-png','-pdf','-m2')


