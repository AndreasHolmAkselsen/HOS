function df = dfRankine(z,zc,A,R,mirror)
% Rankin-type vortex of strength A, centre zx and inner radius R.
% % example:
% R = .5;
% A = 1i;
% zc = .4-.8i;
% z = (-1:.05:1) + 1i*(-2:.05:2).' + zc;
% df = current.dfRankine(z,zc,A,R);
% figure; quiver(real(z),imag(z),real(df),-imag(df))
% hold on; axis equal; plot(R*cos(0:.01:2*pi)+real(zc),R*sin(0:.01:2*pi)+imag(zc),'k')
% % or, with two vortices
% A = [1i,-1i];
% zc= [.4+.3i,0+.3i];
% df = current.dfRankine(z,zc,A,R); figure; quiver(real(z),imag(z),real(df),-imag(df))

if isvector(A),  A  = shiftdim(A,-2);  end
if isvector(zc), zc = shiftdim(zc,-2); end

r = abs(z-zc);
if nargin < 5 || ~mirror
    df = sum(  A./(z-zc).*(r>=R) + A.*conj(z-zc)./R^2 .*(r<R),3:4);
else
    cjzc = conj(zc);
    cjr = abs(z-cjzc);
    df = sum( A./(z-zc).*(r>=R) -A./(z-cjzc).*(cjr>=R) + A./R^2.*conj((z-zc).*(r<R)-(z-cjzc).*(cjr<R)),3:4);
end