function df = dfRankin(z,zc,A,R)
% Rankin-type vortex of strength A, centre zx and inner radius R.
% % example:
% R = .5;
% A = 1i;
% zc = .4+.3i;
% z = (-1:.05:1) + 1i*(-1:.05:1).' + zc;
% df = dfRankin(z,zc,A,R);
% figure; quiver(real(z),imag(z),real(df),-imag(df))
% hold on; axis equal; plot(R*cos(0:.01:2*pi)+real(zc),R*sin(0:.01:2*pi)+imag(zc),'k')
% % or, with two vortices
% A = [1i,-1i];
% zc= [.4+.3i,0+.3i];
% df = dfRankin(z,zc,A,R); figure; quiver(real(z),imag(z),real(df),-imag(df))

if isvector(A)
    A = shiftdim(A,-2);
    zc = shiftdim(zc,-2);
end
r = abs(z-zc);
df = sum(  A./(z-zc).*(r>=R) + A.*conj(z-zc)./R^2 .*(r<R),3:4);


