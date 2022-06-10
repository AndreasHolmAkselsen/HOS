function [PhiAdd,flapTime] = callSFoWavemakerFunctions(waveMaker,N,Lx,DO_PADDING)
% Linear wavemaker init, copied from SFo git hosm-nwt2d





x=Lx/N*(0:N);
Nz0=waveMaker.Nz;
extZDomainRatio=waveMaker.extZDomainRatio; %Size of the additional vertical domain (divided by h)
Nz=(1+extZDomainRatio)*Nz0;
%     Th=tanh(kz*Lx);
hF=waveMaker.hingeDepth;                     %Depth of flap's rotation point [m]

%Initialize wave maker signal
switch waveMaker.signal{1}.type
    case 'harmonicRamped'    % regular
        dt=waveMaker.signal{1}.dt;                   %Time step to define flap motion
        h = waveMaker.h;
        [flapThetaDeg,flapTime]=BM.initializeFlapAngle_HarmonicRamp(waveMaker);
    case 'specFile'

        waveDef = timsas_r1.data.WaveDefinition;
        waveDef.readFromSpecFile(waveMaker.signal{1}.specFile);
        h = waveDef.waterDepth;
%         assert(h==waveDef.waterDepth,'Inconsistant water depths. Map: h=%g, spec file: h=%g.',h,waveDef.waterDepth);
        dt = waveDef.dt;                   %Time step to define flap motion
        ts = waveDef.getTimeRealization([0;0]); % returns a time series object
%         ts.simplifyToConstantTimeStep();%time step is constant
%         ts.trim([0,tEndFull]);


%         % regular wave test
%         Ttest = 1; aTest = .0248; ts.value = aTest*sin((2*pi/Ttest).*ts.getTime);

        ts.taper(waveMaker.signal{1}.tRamp*[1,1]);
        [f,hat_x] = FFTfreq(waveDef.dt,ts.value,'positive');
        k = level1.wave.computeWaveNumber(2*pi*f,h,0,0);
        kh = k*h; th = tanh(kh); ch = cosh(kh);
        wbl = waveMaker.hingeDepth;
        d = max(0,h-wbl);
        cosCosTerm = exp(k*(d-h)).* (1+exp(-2*k*d))./(1+exp(-2*k*h)); % =cosh(k*d)./cosh(k*h) 
%         cosCosTerm = exp(-k*wbl).*(1+exp(-2*k*(h-wbl)))./(1+exp(-2*k*h));
        TF = 2*th./(th+kh./ch.^2).*(th - (1-cosCosTerm)./(k*wbl));
        FFTthetaDeg = hat_x./TF / wbl *180/pi; % to linear order!
        FFTthetaDeg(1) = 0;
        zero1 = zeros(mod(waveDef.nt-1,2));
        flapThetaDeg = ifft([FFTthetaDeg(1:end);zero1;conj(FFTthetaDeg(end:-1:2))]);
        assert(isreal(flapThetaDeg))
%         flapThetaDeg = [flapThetaDeg;zeros(ceil(waveMaker.signal{1}.tFinalStill/dt),1)];
%         flapTime = (0:length(flapTheta)-1)*dt;
        flapTime = ts.getTime;
%         figure();plot(flapTime,flapThetaDeg,'-');grid on; ylabel('\theta');xlabel('Time [s]')
        
end

[flapX,flapZ]=BM.initializeFlapMotion_SingleFlap(flapThetaDeg,h,hF,Nz0);
[X,z]=BM.extendToAdditionalDomain(flapX,flapZ,extZDomainRatio);

% flapMotion=BM.computeFlapMotionDerivatives(X,dt,z(2)-z(1));
flapMotion=zeros(size(X,1),size(X,2),3); % only d_t and d_tt needed at first order.
[~,flapMotion(:,:,2)]=gradient(X,dt);
[~,flapMotion(:,:,3)]=gradient(flapMotion(:,:,2),dt);



%Precompute the chi-matrix used to obtain the additional WM-potential on z=0

%Fine x-grid for anti-aliaisng
% AHA comment: we have to interpolate anyway, so doing the padding here isn't really necessary, but we'll keep it.
p= 1+3*DO_PADDING; % p=4;
Nd=N*(p+1)/2;
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

PhiAdd=zeros(size(dphi1Add_dx_Z0,2),size(dphi1Add_dx_Z0,1),3);
PhiAdd(:,:,1)=dphi1Add_dx_Z0.';
PhiAdd(:,:,2)=dphi1Add_dz_Z0.';
PhiAdd(:,:,3)=dphi1Add_dt_Z0.';
flapTime = flapTime.';

% PhiAdd=zeros(size(dphi1Add_dx_Z0,1),size(dphi1Add_dx_Z0,2),3);
% PhiAdd(:,:,1)=dphi1Add_dx_Z0;
% PhiAdd(:,:,2)=dphi1Add_dz_Z0;
% PhiAdd(:,:,3)=dphi1Add_dt_Z0;
% 
% 
% % AHA: flip!
% PhiAdd = permute(PhiAdd,[2,1,3]); 
% flapTime = flapTime.';
