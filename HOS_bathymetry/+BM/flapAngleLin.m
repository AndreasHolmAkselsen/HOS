function theta = flapAngleLin(T,ka,h,wbl)

l=wbl-h;
d = min(h,max(-l,0));
k = findWaveNumbers(2*pi/T,h,0,0);

Lambda = .5*( k*h./cosh(k*h).^2 + tanh(k*h) );
Gamma =  tanh(k*h) - 1./(k*(h+l)).*(1-cosh(k*d)./cosh(k*h));
c =  tanh(k*h).*Gamma./Lambda;
% A = ka/k;
% hX = A./(1i.*c);
% theta = abs(hX/wbl);

theta =  ka/(k*c*wbl) * 180/pi;