function [z,df,zz] = fz_volumeConverving(xx_,yy_,theta,h,wbl,wbOverWater)
    assert(all(diff(xx_)>0)&&all(diff(yy_)>0),'increasing zz vectors assumed');
    xx = fliplr(xx_); % integrating from right-to-left
    yy = yy_+h; 
    
    d = wbl+wbOverWater;
    H = h+wbl+2*wbOverWater;
    [z,df,zz] = fz0(xx,yy,theta,H,d);
    
    z  = fliplr(z ) - 1i*h - real(z(1,end,:));
    df = fliplr(df);
    zz = fliplr(zz)-1i*h;
end
function [z,df,zz] = fz0(xx,yy,theta,H,d)
zz = xx + 1i*yy;
sig0 = H + [-2;-1;1;2]*d;

sig = zeros(4,1,numel(theta));
for i = 1:length(theta)
    sig(:,:,i) = LMFnlsq(@(sig) fixHingeDepth(sig,theta(i),H,d),sig0); % {LMFnlsq,newtonraphson,fminsearch}
end

zz2 = zz.^2;
df = ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);
df_yh = .5*(df(1:end-1,1,:)+df(2:end,1,:));
z1 = cumsum([zeros(size(theta));df_yh.*diff(1i*yy)],1);
df_xh = .5*(df(:,1:end-1,:)+df(:,2:end,:));
z =  cumsum([z1,df_xh.*diff(xx)],2) ;
end
function err = fixHingeDepth(sig,theta,H,d)
d_xi = 1e-6;
df_zz2 = @(zz2) ((zz2+sig(1)^2)./(zz2+sig(4)^2).*((zz2+sig(3)^2)./(zz2+sig(2)^2)).^2).^(theta/pi);
dy = @(yy) real( df_zz2((1i*yy+d_xi).^2) ); % real beacuse it is integrated by dzz = 1i*dyy;
y = zeros(5,1); sig0 = [0;sig];
for i = 1:4
    y(i+1) = y(i) + integral( dy,sig0(i),sig0(i+1));
end
err = y(2:5) - (H+[-2;-1;1;2]*d);
end