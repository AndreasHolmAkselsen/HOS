clear variables
close all
cCol={'b','r','g','k','m'};
simPath='.\..\simResults';

%HOS order
M=3:6;

%Main loop
for iM=1:length(M)
    
    %Read simulation data
    specFile=['testMELHOS_RegWaveM',num2str(M(iM)),'_eps0.1'];   
    load([simPath '\' specFile,'_velNWT.mat'],'-MAT');
    
    %Extract velocity profile under wave crest
    if iM==1
        uProf=zeros(length(zOut),length(M));
    end
    for iz=1:length(uProf)
        uProf(iz,iM)=max(max(u(:,iz,:)));
    end
end

%Compute 5th order velocity profile from Fenton(1985)
load([simPath '\' specFile,'_simNWT.mat'],'-MAT');
k0=simRes.nwtSpec.init.waveComp{1}.nk0*2*pi/simRes.nwtSpec.sim.Lx;
sig0=sqrt(k0*9.81);
ep=simRes.nwtSpec.init.waveComp{1}.a*k0;
epE=ep*(1+ep^2/2-ep^4/2);
etamax=(epE+epE^2/2+2*epE^4/3)/k0;
uProfFenton=sig0/k0*((epE-epE^3/2)*exp(k0*zOut)+epE^4*exp(2*k0*zOut));
% uProfFenton=sig0/k0*((epE-epE^3/2)*exp(k0*zOut)+epE^4*exp(2*k0*zOut)+epE^5/24*(-37*exp(k0*zOut)+6*exp(3*k0*zOut)));
% etamax=(epE+epE^2/2+2)/k0;
% uProfFenton=sig0/k0*((epE-epE^3/2)*exp(k0*zOut));
uProfFenton(zOut>etamax)=NaN;

figure
hold on
plot(uProfFenton*k0/sig0,zOut*k0,cCol{1},'linewidth',2)
for iM=1:length(M)
    plot(uProf(:,iM)*k0/sig0,zOut*k0,cCol{iM+1},'linewidth',1)
end
grid on
xlim([0.08 0.115])
ylim([-0.2 0.11])