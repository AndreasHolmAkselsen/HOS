function dYdt=dXZPhiLS_dt_ExpRampO5(t,Y,k,M,h,HLP,Tramp,nRamp)

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
[U,W,~]=computeUW_detailed(deltaX,deltaZ,AphiS,k,M,h);
if size(U,1)<=5
    U=[U;zeros(6-size(U,1),size(U,2))];
    W=[W;zeros(6-size(W,1),size(W,2))];
end

%Anti-aliasing treatment
U_AA=zeroPadding(U,Nd);
W_AA=zeroPadding(W,Nd);

%Compute squared velocities up to order 5
U2_15_AA=zeros(1,size(U_AA,2));
W2_15_AA=zeros(1,size(W_AA,2));
for m=2:5
    for i=1:m-1
        U2_15_AA=U2_15_AA+U_AA(i,:).*U_AA(m-i,:);
        W2_15_AA=W2_15_AA+W_AA(i,:).*W_AA(m-i,:);
    end
end

%Compute time-derivative of potential at particle position with anti-aliasing treatment
dphiLS_dt_AA=0.5*(U2_15_AA+W2_15_AA+F*(sum(U_AA,1).^2+sum(W_AA,1).^2-U2_15_AA-W2_15_AA));
AdphiLS_dt_AA=getModeAmplitudes(dphiLS_dt_AA);
% dphiLS_dt=N2/2/Nd*getFunctionFromModeAmplitudes([0,AdphiLS_dt_AA(2:N2/2)])-g*deltaZ+mean(dphiLS_dt_AA);
dphiLS_dt=N2/2/Nd*getFunctionFromModeAmplitudes([0,AdphiLS_dt_AA(2:N2/2)])-g*deltaZ;

%Low pass filter the time derivative to remove HF-numerical noise
Adxs_dt=HLP.*getModeAmplitudes(sum(U(1:5,:),1)+F*sum(U(6:end,:),1));
Adzs_dt=HLP.*getModeAmplitudes(sum(W(1:5,:),1)+F*sum(W(6:end,:),1));
AdphiLS_dt=HLP.*getModeAmplitudes(dphiLS_dt);
% Adxs_dt(abs(Adxs_dt)<1e-10)=0;
% Adzs_dt(abs(Adzs_dt)<1e-10)=0;
% AdphiLS_dt(abs(AdphiLS_dt)<1e-10)=0;
dxs_dt=getFunctionFromModeAmplitudes(Adxs_dt);
dzs_dt=getFunctionFromModeAmplitudes(Adzs_dt);
dphiLS_dt=getFunctionFromModeAmplitudes(AdphiLS_dt);

%Return time derivative for ode-solver
dYdt=[dxs_dt,dzs_dt,dphiLS_dt].';

% if t>2
%     keyboard
% end

end
