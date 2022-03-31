function param = setCurrentParam(param)
% moved to separat function for pre-computation.
h = param.basin.depth;
d = param.current.wallDepth;
U0 = param.current.U0;
dk = pi/h;
% k = shiftdim(1:param.current.nk,-1)*dk;
k = (1:param.current.nk)*dk;
Lambda = .5*k*h + .25*sin(2*k*h);
switch param.current.inletType
    case 'uniform'
        UI =  U0/(1-d/h);
        Gamma = - UI.*sin(k*d)./k;
    case 'linear'
        UI = 2*U0./(1-d/h);
        Gamma =  UI.*(cos(k*d)-cos(k*h)-(h-d)*k.*sin(k*d))./(k.^2*(h-d));
    case 'flipLinear'
        UI = 2*U0./(1-d/h);
        Gamma = UI.*(cos(k*h)-cos(k*d))./(k.^2.*(h-d));
end
param.current.k = shiftdim(k,-1);
param.current.hphi = shiftdim(-Gamma./Lambda,-1);
end