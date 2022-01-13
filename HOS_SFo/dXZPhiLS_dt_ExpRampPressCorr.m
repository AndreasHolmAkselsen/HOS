function dYdt=dXZPhiLS_dt_ExpRampPressCorr(t,Y,k,M,h,HLP,Tramp,nRamp,abcd)

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

%Compute z-correction to fullfil mass conservation + p=cst on free surface
AphiS=getModeAmplitudes(phiLS);
deltaZCorr=zS+computeZPressCorr(deltaX,zS,AphiS,k,M,h,abcd,F);

%Recompute particle velocities with corrected particle height
[UlinCorr,UnlCorr,WlinCorr,WnlCorr,~]=computeUW(deltaX,deltaZCorr,AphiS,k,M,h);

%Recompute time-derivative of potential with corrected particle height
UCorr_AA=zeroPadding(UlinCorr+UnlCorr,Nd);
WCorr_AA=zeroPadding(WlinCorr+WnlCorr,Nd);
U2Corr_AA=(UCorr_AA.^2+WCorr_AA.^2)*F;
AU2Corr_AA=getModeAmplitudes(U2Corr_AA);
U2Corr=N2/2/Nd*getFunctionFromModeAmplitudes([0,AU2Corr_AA(2:N2/2)]);
DphiLSCorr_dt=0.5*U2Corr-g*deltaZCorr+0.5*mean(U2Corr_AA);

%Low pass filter the time derivative remove HF-numerical noise
dxs_dt=getFunctionFromModeAmplitudes(HLP.*getModeAmplitudes(UlinCorr+F*UnlCorr));
dzs_dt=getFunctionFromModeAmplitudes(HLP.*getModeAmplitudes(WlinCorr+F*WnlCorr));
DphiLSCorr_dt=getFunctionFromModeAmplitudes(HLP.*getModeAmplitudes(DphiLSCorr_dt));

%Return time derivative for ode-solver
dYdt=[dxs_dt,dzs_dt,DphiLSCorr_dt].';

end
