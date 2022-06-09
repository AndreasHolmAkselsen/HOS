function [theta,t]=initializeFlapAngle_HarmonicRamp(waveMaker)

dt=waveMaker.signal{1}.dt;                   %Time step to define flap motion
thetaAmp=waveMaker.signal{1}.thetaAmpDeg;    %Amplitude of the single flap motion [deg]
T=waveMaker.signal{1}.T;                     %Period of the single flap motion [s]
tMax=waveMaker.signal{1}.tEnd;               %Duration of flap motion signal
tRamp=waveMaker.signal{1}.tRamp;             %Duration of ramp (beginning and end) [s]


t=(0:dt:tMax-1).';

%Compute flap angle without ramp
theta=thetaAmp*sin(2*pi*t/T);

%Add cosine-ramps at the beginning and at the end of the signal
tapWin=ones(size(theta));
tapWin(t<tRamp)=0.5*(1-cos(pi*t(t<tRamp)/tRamp));
tapWin(t>(tMax-tRamp))=0.5*(1+cos(pi*(t(t>(tMax-tRamp))-(tMax-tRamp))/tRamp));
theta=theta.*tapWin;


end





% %Compute time vector
% t=(0:dt:(tMax+tFinalStill)).';
% tFlap=(0:dt:tMax).';
% 
% %Compute flap angle without ramp
% theta=thetaAmp*sin(2*pi*tFlap/T);
% 
% %Add cosine-ramps at the beginning and at the end of the signal
% tapWin=ones(size(theta));
% tapWin(tFlap<tRamp)=0.5*(1-cos(pi*tFlap(tFlap<tRamp)/tRamp));
% tapWin(tFlap>(tMax-tRamp))=0.5*(1+cos(pi*(tFlap(tFlap>(tMax-tRamp))-(tMax-tRamp))/tRamp));
% theta=theta.*tapWin;
% 
% %Add still signal at the end
% theta=[theta;zeros(length(t)-length(tFlap),1)];