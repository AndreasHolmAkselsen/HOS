function Y_t = HOS_Taylor_closed(t,Y,param)

global timeReached
[vphiS,eta] = deal(Y(1:end/2,:),Y(end/2+1:end,:));

if t-timeReached > 1
    timeReached = floor(t);
    fprintf('Time %ds reached (%.1f%%)\n',timeReached,100*timeReached/param.t_end);
end

wNl = param.nonLinRamp(t);

N = size(eta,1);
dxi = param.map.xi(2)-param.map.xi(1);
dk = pi/(dxi*(N-1));
kx = (0:N-1)'*dk;
k = kx;

% w is the vertical velocity in the zeta-plane; \varphi_\sigma
[w_lin,w_nl] = phiComponentsHOS_closed(vphiS,eta,k,param);

FFTeta = cosfft(eta);
FFTvphiS = cosfft(vphiS);
try
    zzS  = param.map.xi+1i*eta;
    h    = param.map.fy(zzS);
    f_zz = param.map.f_zz(zzS);
catch ME
    warning(ME.identifier,'%s',ME.message);
    Y_t = nan(size(Y));
    return
end
JInv = abs(f_zz).^(-2);

if param.DO_PADDING
    error('padding not yet supported in closed domain simulation. Check whether it is straight forward.')
    Nd = N*(4+1)/2;
    w_lin = ifft(fftPad(fft(w_lin),Nd));
    w_nl = ifft(fftPad(fft(w_nl),Nd));
    h = ifft(fftPad(fft(h),Nd));
    JInv = ifft(fftPad(fft(JInv),Nd));
else
    Nd = N;
end
    
eta_xi  =  isinfft(fftPad(-kx.*FFTeta,Nd));
vphiS_xi =  isinfft(fftPad(-kx.*FFTvphiS,Nd));
w = w_lin+w_nl;

eta_t  =  JInv.*( w_lin + wNl.*(  w_nl  + eta_xi.^2.*w - vphiS_xi.*eta_xi ) );
vphiS_t = - h*param.g + JInv.* wNl.*( -.5*vphiS_xi.^2  + .5*(1+eta_xi.^2).*w.^2 );

if isfield(param,'beach')
%     phi_xi = vphiS_xi-w.*eta_xi;
%     U = JInv.*((phi_xi+1i*w).*f_zz);
%     h_x = eta_xi - imag(f_zz)./real(f_zz);
%     vphiS_t = vphiS_t - param.beach.nu.*( imag(U)-h_x.*real(U) );
%     figure, plot(param.map.xi,imag(U)-h_x.*real(U),param.map.xi,real(f_zz).*eta_t,'--')
    vphiS_t = vphiS_t - param.beach.nu.* real(f_zz).*eta_t; % approximate version
end
% Add waveMaker, if any (copy from SFo git hosm-nwt2d)
if isfield(param,'waveMaker') && ~(t<param.waveMaker.time(1)||t>=param.waveMaker.time(end-1))
    %Interpolate partial derivatives of additional potential
    % param.waveMaker.phiAdd(i,n,j) are: j-{x,y,t} derivatives; i: xi_i positions; n: times n.
    dtFlap=param.waveMaker.time(2)-param.waveMaker.time(1);
    indTime=1+fix((t-param.waveMaker.time(1))/dtFlap);
    w1=(param.waveMaker.time(indTime+1)-t)/dtFlap;
    w2=(t-param.waveMaker.time(indTime))/dtFlap;
    phiAdd_x=w1*param.waveMaker.phiAdd(:,indTime,1)+w2*param.waveMaker.phiAdd(:,indTime+1,1);
    phiAdd_y=w1*param.waveMaker.phiAdd(:,indTime,2)+w2*param.waveMaker.phiAdd(:,indTime+1,2);
    phiAdd_t=w1*param.waveMaker.phiAdd(:,indTime,3)+w2*param.waveMaker.phiAdd(:,indTime+1,3);
    phiAdd_y=phiAdd_y-mean(phiAdd_y);
    
    Uadd__f_zz = (phiAdd_x+1i*phiAdd_y)./f_zz;
    eta_t = eta_t - eta_xi.*real(Uadd__f_zz) + imag(Uadd__f_zz) ;
    vphiS_t = vphiS_t - phiAdd_t - vphiS_xi.*real(Uadd__f_zz) -0.5*(phiAdd_x.^2+phiAdd_y.^2);    
end

% Unpad, lowpass filter and dampen:
M = ceil(N/2)-1;
M = N;
Md = param.kd__kmax*M;
mu = param.rDampingDim*M*(((0:N-1)'-Md)/(M-Md)).^2.*((0:N-1)'>Md);    
% kMax = (N-1)*dk;
% kd = param.kd__kmax*kMax*dk;
% mu = param.rDamping*kMax*((k-kd)/(kMax-kd)).^2.*(k>kd);     
kFilter = k<=param.iModeCut*dk;  
Y_t = real([ icosfft(kFilter.*fftPad(cosfft(vphiS_t),N)-mu.*FFTvphiS)
             icosfft(kFilter.*fftPad(cosfft(eta_t  ),N)-mu.*FFTeta  ) ]);        
end

