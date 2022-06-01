function [z,f_zz,zz] = fz_volumeConverving(xx_,yy_,theta,h,wbl,wbOverWater)
    assert(all(diff(xx_)>0)&&all(diff(yy_)>0),'increasing zz vectors assumed');
    xx = fliplr(xx_); % integrating from right-to-left
    yy = yy_+h; 
    
    d = wbl+wbOverWater;
    H = h+wbl+2*wbOverWater;
    [z,f_zz,zz] = fz0(xx,yy,theta,H,d);
    
    z  = fliplr(z ) - 1i*h - real(z(1,end,:));
    f_zz = fliplr(f_zz);
    zz = fliplr(zz)-1i*h;
end
function [z,f_zz,zz] = fz0(xx,yy,theta,H,d)
zz = xx + 1i*yy;
sig0 = H + [-2;-1;1;2]*d;

sig = zeros(4,1,numel(theta));
for i = 1:length(theta)
    sig(:,:,i) = LMFnlsq(@(sig) fixHingeDepth(sig,theta(i),H,d),sig0); % {LMFnlsq,newtonraphson,fminsearch}
end

zz2 = zz.^2;
f_zz = ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);
f_zz_yh = .5*(f_zz(1:end-1,1,:)+f_zz(2:end,1,:));
z1 = cumsum([zeros(size(theta));f_zz_yh.*diff(1i*yy)],1);
f_zz_xh = .5*(f_zz(:,1:end-1,:)+f_zz(:,2:end,:));
z =  cumsum([z1,f_zz_xh.*diff(xx)],2) ;
end
function err = fixHingeDepth(sig,theta,H,d)
d_xi = 1e-6;
f_zz_zz2 = @(zz2) ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);
dy = @(yy) real( f_zz_zz2((1i*yy+d_xi).^2) ); % real beacuse it is integrated by dzz = 1i*dyy;
y = zeros(5,1); sig0 = [0;sig];
for i = 1:4
    y(i+1) = y(i) + integral( dy,sig0(i),sig0(i+1));
end
err = y(2:5) - (H+[-2;-1;1;2]*d);
end