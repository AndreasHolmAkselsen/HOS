function [PhiAdd,flapTime] = callSFoWavemakerFunctions(waveMaker,Nx,Lx,h,DO_PADDING)
% Linear wavemaker init, copied from SFo git hosm-nwt2d

x=Lx/Nx*(0:Nx);
Nz0=waveMaker.Nz;
extZDomainRatio=waveMaker.extZDomainRatio; %Size of the additional vertical domain (divided by h)
Nz=(1+extZDomainRatio)*Nz0;
%     Th=tanh(kz*Lx);
hF=waveMaker.hingeDepth;                     %Depth of flap's rotation point [m]
thetaAmp=waveMaker.signal{1}.thetaAmpDeg;    %Amplitude of the single flap motion [deg]
T=waveMaker.signal{1}.T;                     %Period of the single flap motion [s]
dt=waveMaker.signal{1}.dt;                   %Time step to define flap motion
tMax=waveMaker.signal{1}.tEnd;               %Duration of flap motion signal
tRamp=waveMaker.signal{1}.tRamp;             %Duration of ramp (beginning and end) [s]
tFinalStill=waveMaker.signal{1}.tFinalStill; %Duration of the still period at the end
%Initialize wave maker signal
[flapTheta,flapTime]=BM.initializeFlapAngle_HarmonicRamp(thetaAmp,dt,tMax,T,tRamp,tFinalStill);
[flapX,flapZ]=BM.initializeFlapMotion_SingleFlap(flapTheta,h,hF,Nz0);
[X,z]=BM.extendToAdditionalDomain(flapX,flapZ,extZDomainRatio);
flapMotion=BM.computeFlapMotionDerivatives(X,dt,z(2)-z(1));



%Precompute the chi-matrix used to obtain the additional WM-potential on z=0
%Fine x-grid for anti-aliaisng
p= 1+3*DO_PADDING; % p=4;
Nd=Nx*(p+1)/2;
x_AA=(0:Nd)/Nd*Lx;


kz=pi/(h*(1+extZDomainRatio))*(0:Nz);
C=kz.'*(Lx-x_AA);
D=kz.'*repmat(Lx,1,length(x_AA));
limMask=C<10;
chi=exp(C-D);   %Simplified expression to avoid infty/infty terms
chi(limMask)=cosh(C(limMask))./cosh(D(limMask));
shi=chi;
shi(limMask)=sinh(C(limMask))./cosh(D(limMask));
[xMat,kMat]=meshgrid(x_AA,kz);
xChiMin=min(xMat((chi<=10^-4)&(xMat/Lx)>(kMat/max(kz))));
kChiMin=min(kMat((chi<=10^-4)&(xMat/Lx)<(kMat/max(kz))));
chiLargeX=chi(kz<=kChiMin,:);
chiLargeK=chi(kz>kChiMin,x<=xChiMin);
xShiMin=min(xMat((shi<=10^-4)&(xMat/Lx)>(kMat/max(kz))));
kShiMin=min(kMat((shi<=10^-4)&(xMat/Lx)<(kMat/max(kz))));
shiLargeX=shi(kz<=kShiMin,:);
shiLargeK=shi(kz>kShiMin,x<=xShiMin);


%Pre-compute phi1Add from flap velocity
Th=tanh(kz*Lx);
dphi1Add_dx=flapMotion(:,:,2);
%     Adphi1Add_dx=getCosModeAmplitudes(dphi1Add_dx);
Adphi1Add_dx=cosfft(dphi1Add_dx.').';
Aphi1Add=BM.getXAddDerivativeX0Amp(Adphi1Add_dx,-1,kz,Th);
Adphi1Add_dz=BM.getZAddDerivativeX0Amp(Aphi1Add,1,kz);
dphi1Add_dxdt=flapMotion(:,:,3);
%     Adphi1Add_dxdt=getCosModeAmplitudes(dphi1Add_dxdt);
Adphi1Add_dxdt=cosfft(dphi1Add_dxdt.').';
Adphi1Add_dt=BM.getXAddDerivativeX0Amp(Adphi1Add_dxdt,-1,kz,Th);

%Pre-compute partial derivatives of phi1Add at z=0
dphi1Add_dx_Z0=zeros(size(Adphi1Add_dx,1),length(x_AA));
dphi1Add_dz_Z0=zeros(size(Adphi1Add_dx,1),length(x_AA));
dphi1Add_dt_Z0=zeros(size(Adphi1Add_dx,1),length(x_AA));
%         tic
for iT=1:size(Adphi1Add_dx,1)
    indKeep=abs(Adphi1Add_dx(iT,:))>1e-4*max(abs(Adphi1Add_dx(iT,:)));
    dphi1Add_dx_Z0(iT,:)=BM.getAddFunctionAtZ0_cos(Adphi1Add_dx(iT,:),kz,h,shiLargeX,shiLargeK,indKeep);
    indKeep=abs(Adphi1Add_dz(iT,:))>1e-4*max(abs(Adphi1Add_dz(iT,:)));
    dphi1Add_dz_Z0(iT,:)=BM.getAddFunctionAtZ0_sin(Adphi1Add_dz(iT,:),kz,h,chiLargeX,chiLargeK,indKeep);
    indKeep=abs(Adphi1Add_dt(iT,:))>1e-4*max(abs(Adphi1Add_dt(iT,:)));
    dphi1Add_dt_Z0(iT,:)=BM.getAddFunctionAtZ0_cos(Adphi1Add_dt(iT,:),kz,h,chiLargeX,chiLargeK,indKeep);
end
%         toc
PhiAdd=zeros(size(dphi1Add_dx_Z0,1),size(dphi1Add_dx_Z0,2),3);
PhiAdd(:,:,1)=dphi1Add_dx_Z0;
PhiAdd(:,:,2)=dphi1Add_dz_Z0;
PhiAdd(:,:,3)=dphi1Add_dt_Z0;


% AHA: flip!
PhiAdd = permute(PhiAdd,[2,1,3]); 
flapTime = flapTime.';
