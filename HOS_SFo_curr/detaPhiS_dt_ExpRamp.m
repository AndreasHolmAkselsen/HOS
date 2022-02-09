function dYdt=detaPhiS_dt_ExpRamp(t,Y,k,M,h,HLP,Tramp,nRamp,df,deAliasingProductOrder)

%Extract surface elevation and potential
g=9.81;
N2=length(Y)/2;
eta=Y(1:N2).';
phiS=Y(N2+1:end).';

%Compute ramp function
if t<0
    F=0;
else
    F=1-exp(-(t/Tramp)^nRamp);
end

%Compute Fourier amplitudes of eta and phiS
Aeta=getModeAmplitudes(eta);
AphiS=getModeAmplitudes(phiS);

%Find extended number of modes for anti-aliasing (4th order product)
N=length(Aeta);
% p=4;
p = deAliasingProductOrder;
Nd=N*(p+1)/2;

%Compute space gradients with extended number of modes + low pass filter 
Aeta_AA=Nd/N*[Aeta,zeros(1,Nd-N)];
AphiS_AA=Nd/N*[AphiS,zeros(1,Nd-N)];
k_AA=(k(2)-k(1))*(0:(Nd-1));
dAeta_AA=1i*k_AA.*Aeta_AA;
dAphiS_AA=1i*k_AA.*AphiS_AA;
eta_AA=ifft([Aeta_AA,0,fliplr(conj(Aeta_AA(2:end)))]);
grad_eta_AA=ifft([dAeta_AA,0,fliplr(conj(dAeta_AA(2:end)))]);
grad_phiS_AA=ifft([dAphiS_AA,0,fliplr(conj(dAphiS_AA(2:end)))]);

%Compute vertical velocity and extend mode numbers with zero padding
% if t>130
%     figure;plot(log10(1e-10+abs(Aeta)));grid on
%     keyboard
% end
[W,~]=computeW_detailed(eta,AphiS,k,M,h);
W_lin_AA=zeroPadding(W(1,:),Nd);
W_nl_AA=zeroPadding(sum(W(2:end,:),1),Nd);
W_tot_AA=W_lin_AA+W_nl_AA;

% Potential current function
FCurr = F;
L = 2*pi/k(2);

dx_AA = L/(2*Nd-1);
x_AA = (0:2*Nd-1)*dx_AA;
dfS_AA = df(x_AA + 1i*eta_AA);
Phi_x_AA =  real(dfS_AA);
Phi_z_AA = -imag(dfS_AA);

% dx = L/(N2-1);
% x = (0:N2-1)*dx;
% dfS = df(x + 1i*eta);
% Phi_x_AA = zeroPadding( real(dfS),Nd);
% Phi_z_AA = zeroPadding(-imag(dfS),Nd);

%Compute time derivatives and keep only N modes (full anti-aliasing treatment)  
Adeta_dt_AA  = getModeAmplitudes(W_lin_AA + F* (  W_nl_AA     + grad_eta_AA.^2.*W_tot_AA - grad_eta_AA.*grad_phiS_AA) - FCurr.*(Phi_x_AA.*grad_eta_AA - Phi_z_AA)   );
AdphiS_dt_AA = getModeAmplitudes(-g*eta_AA - 0.5*F* grad_phiS_AA.^2 + 0.5*F* (1+grad_eta_AA.^2).*W_tot_AA.^2          - FCurr.*( Phi_x_AA.*grad_phiS_AA + .5*Phi_x_AA.^2 + .5*Phi_z_AA.^2)       );    

%Low pass filter the time derivative to remove HF-numerical noise
Adeta_dt=N/Nd*HLP.*[0,Adeta_dt_AA(2:N)];
AdphiS_dt=N/Nd*HLP.*[0,AdphiS_dt_AA(2:N)];
Adeta_dt(abs(Adeta_dt)<1e-6)=0;
AdphiS_dt(abs(AdphiS_dt)<1e-6)=0;
deta_dt=getFunctionFromModeAmplitudes(Adeta_dt);
dphiS_dt=getFunctionFromModeAmplitudes(AdphiS_dt);
% deta_dt=N/Nd*getFunctionFromModeAmplitudes();
% dphiS_dt=N/Nd*getFunctionFromModeAmplitudes(HLP.*[0,AdphiS_dt_AA(2:N)]);

%Return time derivative for ode-solver
dYdt=[deta_dt.';dphiS_dt.'];


end
