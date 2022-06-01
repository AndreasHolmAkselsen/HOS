function [z,f_zz,zz] = fz_yySymmetric(xx,yy,theta,h,wbl,wbOverWater)
assert(all(diff(xx)>0)&&all(diff(yy)>0),'increasing zz vectors assumed');
xx = fliplr(xx); % integrating from right-to-left

thp = theta/pi;
% d = sqrt(pi).*(wbl+wbOverWater).*sec(theta)./(gamma(1+thp)*gamma(.5-thp));
zz = xx + 1i*yy;
% f_zz = (1+d^2./zz.^2).^thp; % k=1

% stretched variables
h_   = h   + wbOverWater;
wbl_ = wbl + wbOverWater;
zz_  = zz  - 1i*wbOverWater;

% d = fixHingeDepth(D/H,theta)*H;
d = 0*theta;
for i = 1:length(theta)
    d(i) = fzero(@(d__h) fixHingeDepth_fzero(d__h,wbl_/h_,theta(i)),wbl_/h_)*h_;
end
% figure, plot(theta(:),d(:)-linspace(d(1),d(end),length(d))','.-')

f_zz = (1+csch(pi*zz_./(2*h_)).^2.*sin(pi*d./(2*h_)).^2).^thp; 
f_zz_yh = .5*(f_zz(1:end-1,1,:)+f_zz(2:end,1,:));
z1 = cumsum([zeros(size(theta));f_zz_yh.*diff(1i*yy)],1);
z1 = z1-z1(yy==0);% + 1i*wbOverWater;
f_zz_xh = .5*(f_zz(:,1:end-1,:)+f_zz(:,2:end,:));
z =  cumsum([z1,f_zz_xh.*diff(xx)],2) ;
z = z-real(z(1,end,:));

z  = fliplr(z );
f_zz = fliplr(f_zz);
zz = fliplr(zz);

end

function err = fixHingeDepth_fzero(d__h,wbl__h,theta)
% there's a singularity at zz=-1i*d if theta < 0
% either stop a yy=-d-delta_singularity
% or shift integration path d_xi to the right (into the domain):
% delta_singularity = 1e-6*(theta<0); d_xi = 0;
delta_singularity = 0; d_xi = 1e-6*(theta<0);

% nYInt = 10000;
% yi = linspace(-1,-d__h-delta_singularity,nYInt)';
% dL = real( (1-csc(pi/2*(yi-1i*d_xi)).^2.*sin(pi/2*d__h).^2).^(theta/pi) );
% L = sum(.5*(dL(1:end-1)+dL(2:end)) .* diff(yi) )+delta_singularity;

L = integral(@(y) real((1-csc(pi/2*(y-1i*d_xi)).^2.*sin(pi/2*d__h).^2).^(theta/pi)),-1,-d__h-delta_singularity)+delta_singularity;
err = wbl__h-(1-L);
end
