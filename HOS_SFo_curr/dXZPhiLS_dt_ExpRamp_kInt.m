function dYdt=dXZPhiLS_dt_ExpRamp_kInt(t,Y,k,M,h,HLP,Tramp,nRamp)

g=9.81;
N=length(Y)/3;
N2=2*N;
p=2;
Nd=N2*(p+1)/4;

%Extract particle positions and potential
AxS=Y(1:N).';
AzS=Y(N+1:2*N).';
AphiLS=Y(2*N+1:3*N).';
deltaX=getFunctionFromModeAmplitudes([0,AxS(2:end)]);
zS=getFunctionFromModeAmplitudes(AzS);
phiLS=getFunctionFromModeAmplitudes(AphiLS);

%Compute ramp function
if t<0
    F=0;
else
    F=1-exp(-(t/Tramp)^nRamp);
end

%Compute deltaZ, with nonlinear correction to obtain zero mean surface elevation 
eta=sum(computeEta(deltaX,zS,k,M),1);
deltaZ=zS-F*mean(eta);

%Compute particle velocities
AphiS=getModeAmplitudes(phiLS);
[Ulin,Unl,Wlin,Wnl]=computeUW(deltaX,deltaZ,AphiS,k,M,h);

%Anti-aliasing treatment
U_AA=zeroPadding(Ulin+Unl,Nd);
W_AA=zeroPadding(Wlin+Wnl,Nd);

%Compute time-derivative of potential at particle position with anti-aliasing treatment
dphiLS_dt_AA=0.5*(U_AA.^2+W_AA.^2)*F;
AdphiLS_dt_AA=getModeAmplitudes(dphiLS_dt_AA);
dphiLS_dt=N/Nd*getFunctionFromModeAmplitudes([0,AdphiLS_dt_AA(2:N)])-g*deltaZ+mean(dphiLS_dt_AA);

%Low pass filter the time derivative remove HF-numerical noise
dAxs_dt=HLP.*getModeAmplitudes(Ulin+F*Unl);
dAzs_dt=HLP.*getModeAmplitudes(Wlin+F*Wnl);
dAphiLS_dt=HLP.*getModeAmplitudes(dphiLS_dt);

%Return time derivative for ode-solver
dYdt=[dAxs_dt,dAzs_dt,dAphiLS_dt].';

end
