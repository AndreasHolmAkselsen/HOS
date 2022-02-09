function dYdt=dXZPhiLS_dt_AccExpRamp(t,Y,k,M,h,HLP,Tramp,nRamp,tau,abcd)

g=9.81;
N2=(length(Y)-1)/4;
p=2;
Nd=N2*(p+1)/4;

%Extract particle positions and potential
xS=Y(1:N2).';
zS=Y(N2+1:2*N2).';
phiLS=Y(2*N2+1:3*N2).';
dphiLS_dt=Y(3*N2+1:4*N2).';
zm=Y(end);

%Compute ramp function
if t<0
    F=0;
else
    F=1-exp(-(t/Tramp)^nRamp);
end

%Compute deltaX
deltaX=xS-mean(xS);

%Compute mean water level at this time step
AdeltaX=getModeAmplitudes(deltaX);
AddeltaXda=AdeltaX.*(1i*k);
ddeltaXda=getFunctionFromModeAmplitudes(AddeltaXda);
etaBar=mean(zS.*(1+ddeltaXda));
dzm_dt=(etaBar-zm)/tau;

%Set mean water level to zero
deltaZ=zS-F*zm;

%Compute particle velocities
AphiS=getModeAmplitudes(phiLS);
[Ulin,Unl,Wlin,Wnl,AphiE]=computeUW(deltaX,deltaZ,AphiS,k,M,h);

% Anti-aliasing treatment
U_AA=zeroPadding(Ulin+Unl,Nd);
W_AA=zeroPadding(Wlin+Wnl,Nd);

%Compute time-derivative of potential at particle position with anti-aliasing treatment
U2_AA=(U_AA.^2+W_AA.^2);
AU2_AA=getModeAmplitudes(U2_AA);
U2=N2/2/Nd*getFunctionFromModeAmplitudes([0,AU2_AA(2:N2/2)]);
AphiSDot=getModeAmplitudes(-0.5*U2-g*deltaZ-0.5*mean(U2_AA));

%Compute the Eulerian time derivative of the surface velocity
[UDotlin,UDotnl,WDotlin,WDotnl,~]=computeUW(deltaX,deltaZ,AphiSDot,k,M,h);

%Compute acceleration's convective term 
[axConv,azConv]=computeAccConv(deltaX,deltaZ,AphiE,k,M,h,abcd);

%Compute total particle acceleration with antialiasing treatment
ax_AA=zeroPadding(UDotlin+UDotnl+axConv,Nd);
az_AA=zeroPadding(WDotlin+WDotnl+azConv,Nd);

%Compute the double Lagrangian time derivative of phiS
dphiSDot_dt_AA=zeroPadding(-g*(Wlin+F*Wnl),Nd);
dphiSDot_dt_AA=dphiSDot_dt_AA+U_AA.*ax_AA+W_AA.*az_AA;
AdphiSDot_dt_AA=getModeAmplitudes(dphiSDot_dt_AA);
dphiSDot_dt=N2/2/Nd*getFunctionFromModeAmplitudes(AdphiSDot_dt_AA(1:N2/2));

%Low pass filter the time derivative to remove HF-numerical noise
dxs_dt=getFunctionFromModeAmplitudes(HLP.*getModeAmplitudes(Ulin+F*Unl));
dzs_dt=getFunctionFromModeAmplitudes(HLP.*getModeAmplitudes(Wlin+F*Wnl));
dphiLS_dt=getFunctionFromModeAmplitudes(HLP.*getModeAmplitudes(dphiLS_dt));
dphiLSDot_dt=getFunctionFromModeAmplitudes(HLP.*getModeAmplitudes(dphiSDot_dt));

%Return time derivative for ode-solver
dYdt=[dxs_dt,dzs_dt,dphiLS_dt,dphiLSDot_dt,dzm_dt].';

end
