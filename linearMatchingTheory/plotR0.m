clear
% see also c:/unix_data/bathymetryStep/Seb_R0
addpath c:/gits/wavemaker

nEv = 200;


% T = 2;
% h_d = 1.0;
% h_s = .6;

h_d = 6;
h_s = [.5,1,2,3,4];
T = (.5:.05:5)';


R0val = 0*(h_s+T);
for iT = 1:length(T)
    for ihs = 1:length(h_s)
        R0val(iT,ihs) = R0(T(iT),h_d,h_s(ihs),nEv);
    end
end

figure('color','w');
plot(T,R0val,'linewidth',1.5);
legend("h_s = "+h_s'+"m",'location','northwest');
grid on;
xlabel('Periods [s]');ylabel('|R_0|')
xlim(T([1,end]));
ylim([0,.5]);