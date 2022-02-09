function zCorr=computeZPressCorr(deltaX,deltaZ,AphiS,k,M,h,abcd,F)

g=9.81;
N2=length(deltaZ);
p=2;
Nd=N2*(p+1)/4;

%Compute deltaZ, with nonlinear correction to obtain zero mean surface elevation
AdeltaX=getModeAmplitudes(deltaX);
AddeltaXda=AdeltaX.*(1i*k);
ddeltaXda=getFunctionFromModeAmplitudes(AddeltaXda);
etaBar=mean(deltaZ.*(1+ddeltaXda));
deltaZ=deltaZ-F*etaBar;

%Compute particle velocities
[Ulin,Unl,Wlin,Wnl,AphiE]=computeUW(deltaX,deltaZ,AphiS,k,M,h);

%Anti-aliasing treatment
U_AA=zeroPadding(Ulin+Unl,Nd);
W_AA=zeroPadding(Wlin+Wnl,Nd);

%Compute time-derivative of potential at particle position with anti-aliasing treatment
U2_AA=(U_AA.^2+W_AA.^2)*F;
AU2_AA=getModeAmplitudes(U2_AA);
U2=N2/2/Nd*getFunctionFromModeAmplitudes([0,AU2_AA(2:N2/2)]);
dphiLS_dt=-0.5*U2-g*deltaZ-0.5*mean(U2_AA);
AphiSDot=getModeAmplitudes(dphiLS_dt);

%Compute time derivative of the horizontal velocity
[UDotlin,UDotnl,WDotlin,WDotnl,~]=computeUW(deltaX,deltaZ,AphiSDot,k,M,h);

%Compute acceleration's convective term 
[axConv,azConv]=computeAccConv(deltaX,deltaZ,AphiE,k,M,h,abcd);

%Compute total particle acceleration
ax=UDotlin+F*(UDotnl+axConv);
az=WDotlin+F*(WDotnl+azConv);

%Compute a-dependent z-correction to have constant pressure on surface
AdeltaZ=getModeAmplitudes(deltaZ);
AddeltaZda=AdeltaZ.*(1i*k);
ddeltaZda=getFunctionFromModeAmplitudes(AddeltaZda);
dpda=ax.*(1+ddeltaXda)+(az+g).*ddeltaZda;
dzCorrda=-dpda./(az+g);
AdzCorrda=getModeAmplitudes(dzCorrda);
AzCorr=AdzCorrda./(1i*k);
AzCorr(1)=0;
zCorr=getFunctionFromModeAmplitudes(AzCorr);
C=-mean((deltaZ+zCorr).*(1+ddeltaXda));
zCorr=zCorr+C-F*etaBar;

end