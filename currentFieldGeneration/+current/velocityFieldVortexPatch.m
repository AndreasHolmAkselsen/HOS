function [u,w,phi,psi] = velocityFieldVortexPatch(x,z,param)
if ~isfield(param.current,'hphi') || isempty(param.current.hphi)
    param = current.setCurrentParam(param);
end

h = param.basin.depth;
% centre = param.vortexPatch.centre;
% strength =  param.vortexPatch.strength;
% width = param.vortexPatch.size(1);
% height = param.vortexPatch.size(2);
% nSize = param.vortexPatch.nSize;

% zzPatch = linspace(0,width,nSize(2))- width/2 + 1i*(linspace(0,height,nSize(1)).'-height/2);
% % zzPatch( imag(zzPatch)<(real(zzPatch)/width-.5)*height+height/2 )=[];
% zzPatch = centre + zzPatch(:).';
% strength = shiftdim(strength,-1); zzPatch = shiftdim(zzPatch,-1);
% A_j = .5*strength.*abs(imag(centre))/numel(zzPatch);

zzPatch = param.vortexPatch.zzPatch;
A_j = param.vortexPatch.A_j;
% A_j = shiftdim(A_j,-1); zzPatch = shiftdim(zzPatch,-1);


if isfield(param.vortexPatch,'R_Rankine') && param.vortexPatch.R_Rankine > 0
    df00 = @(zz) current.dfRankine(zz,zzPatch,A_j,param.vortexPatch.R_Rankine,true);
else
    df00 = @(zz) sum(A_j./(zz-zzPatch) - A_j./(zz-conj(zzPatch)),3:4);
end
if param.vortexPatch.mirror
    df0 = @(zz) df00(zz) - df00(-zz);
    dfVort = @(zz) df0(zz) + df0(zz+2i*h) + df0(zz-2i*h)  + df0(zz+4i*h);
else
    dfVort = @(zz) df00(zz) - df00(-zz);
end
zz = x+1i*z;
k = param.current.k;
df = param.current.U0 + sum( -k.*param.current.hphi.*exp(-k.*zz),3 ) + dfVort(zz);

u = real(df);
w = -imag(df);

if nargout <= 2, return; end

f00  = @(zz) sum(A_j.*log(zz-zzPatch) - A_j.*log(zz-conj(zzPatch)),3:4);
if param.vortexPatch.mirror
    f0 = @(zz) f00(zz)   + f00(-zz);
    fVort = @(zz) f0(zz) + f0(zz+2i*h) + f0(zz-2i*h)  + f0(zz+4i*h);
else
    fVort = @(zz) f00(zz)   + f00(-zz);
end

f = param.current.U0.*zz + sum( param.current.hphi.*exp(-k.*(x+1i*z)),3 ) + fVort(zz);
phi = real(f);
psi = imag(f);
    




% 
% %% temp
% k = shiftdim(param.current.k,-1);
% phi = shiftdim(param.current.hphi,-1);
% v0 = sum( -k.*phi.*exp(-k.*(x+1i*z)),3 );
% 
% dk = pi/h;
% k = param.current.k;
% kk = (1:(param.current.nk*2+1))*dk;
% pphi = [0,-k.*param.current.hphi.*exp(-k.*x),zeros(1,param.current.nk)];
% v1Temp = ifft(pphi)*numel(kk);
% v1 = v1Temp(2:param.current.nk+1).';