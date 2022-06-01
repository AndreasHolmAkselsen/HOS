function [theta,t]=initializeFlapAngle_HarmonicRamp(thetaAmp,dt,tMax,T,tRamp,tFinalStill)

%Compute time vector
t=(0:dt:(tMax+tFinalStill)).';
tFlap=(0:dt:tMax).';

%Compute flap angle without ramp
theta=thetaAmp*sin(2*pi*tFlap/T);

%Add cosine-ramps at the beginning and at the end of the signal
tapWin=ones(size(theta));
tapWin(tFlap<tRamp)=0.5*(1-cos(pi*tFlap(tFlap<tRamp)/tRamp));
tapWin(tFlap>(tMax-tRamp))=0.5*(1+cos(pi*(tFlap(tFlap>(tMax-tRamp))-(tMax-tRamp))/tRamp));
theta=theta.*tapWin;

%Add still signal at the end
theta=[theta;zeros(length(t)-length(tFlap),1)];

end