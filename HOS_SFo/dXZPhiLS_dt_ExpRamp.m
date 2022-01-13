function dYdt=dXZPhiLS_dt_ExpRamp(t,Y,k,M,h,HLP,Tramp,nRamp)

g=9.81;
N2=length(Y)/3;
p=2;
Nd=N2*(p+1)/4;

%Extract particle positions and potential
xS=Y(1:N2).';
zS=Y(N2+1:2*N2).';
phiLS=Y(2*N2+1:3*N2).';

%Compute ramp function
if t<0
    F=0;
else
    F=1-exp(-(t/Tramp)^nRamp);
end

%Compute deltaX
deltaX=xS-mean(xS);

%Compute deltaZ, with nonlinear correction to obtain zero mean surface elevation 
eta=sum(computeEta(deltaX,zS,k,M),1);
deltaZ=zS-F*mean(eta);

%Compute particle velocities
AphiS=getModeAmplitudes(phiLS);
[Ulin,Unl,Wlin,Wnl,~]=computeUW(deltaX,deltaZ,AphiS,k,M,h);

%Anti-aliasing treatment
U_AA=zeroPadding(Ulin+Unl,Nd);
W_AA=zeroPadding(Wlin+Wnl,Nd);

%Compute time-derivative of potential at particle position with anti-aliasing treatment
dphiLS_dt_AA=0.5*(U_AA.^2+W_AA.^2)*F;
AdphiLS_dt_AA=getModeAmplitudes(dphiLS_dt_AA);
dphiLS_dt=N2/2/Nd*getFunctionFromModeAmplitudes([0,AdphiLS_dt_AA(2:N2/2)])-g*deltaZ+mean(dphiLS_dt_AA);

%Low pass filter the time derivative to remove HF-numerical noise
Adxs_dt=HLP.*getModeAmplitudes(Ulin+F*Unl);
Adzs_dt=HLP.*getModeAmplitudes(Wlin+F*Wnl);
AdphiLS_dt=HLP.*getModeAmplitudes(dphiLS_dt);
Adxs_dt(abs(Adxs_dt)<1e-12)=0;
Adzs_dt(abs(Adzs_dt)<1e-12)=0;
AdphiLS_dt(abs(AdphiLS_dt)<1e-12)=0;
dxs_dt=getFunctionFromModeAmplitudes(Adxs_dt);
dzs_dt=getFunctionFromModeAmplitudes(Adzs_dt);
dphiLS_dt=getFunctionFromModeAmplitudes(AdphiLS_dt);

%Return time derivative for ode-solver
dYdt=[dxs_dt,dzs_dt,dphiLS_dt].';

% if t>405
%     keyboard
% end

end
