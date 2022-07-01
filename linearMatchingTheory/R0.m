function val = R0(T,hL,hR,nEv)

% g = 9.81;

hphiI = 1; % irrelevant for value of R0

N = nEv+1;
w = 2*pi/T;
kiR = findWaveNumbers(w,hR,0,nEv);
kiR = abs(real(kiR)) + 1i*(abs(imag(kiR)));
kiL = findWaveNumbers(w,hL,0,nEv);
kiL = abs(real(kiL)) + 1i*(abs(imag(kiL)));

A = nan(2*N);
b = nan(2*N,1);

d_i1 = zeros(N,1);d_i1(1) = 1;

A(1:N,1:N) = eye(N);
A(1:N,N+1:2*N) = kiR.'./kiL.*Gamma(kiL,hL,kiR.',hR)./Lambda(kiL,hL);
b(1:N) = hphiI.*d_i1;

A(N+1:2*N,N+1:2*N) = eye(N);
A(N+1:2*N,1:N) = -Gamma(kiL.',hL,kiR,hR)./Lambda(kiR,hR);

b(N+1:2*N) = hphiI*Gamma(kiL(1),hL,kiR,hR)./Lambda(kiR,hR);

hphi = A\b;
% hphiL = hphi(1:N);
% hphiR = hphi(N+1:2*N);

val = abs(hphi(1)/hphiI);
end

function res = Lambda(k,h)
    res = .5*(h.*sech(k.*h).^2+tanh(k.*h)./k);
end
function res = Gamma(k1,h1,k2,h2)
    res =  (k1.*tanh(k1.*h1)-k2.*tanh(k2.*h2)-k1.*sinh(k1*(h1-h2)).*sech(k1.*h1).*sech(k2.*h2))./(k1.^2-k2.^2);
%     ii = k1==k2;
    ii = abs(k1-k2)<1e-9*(k1+k2);
    if any(ii,'all')
        k = k1+0*k2;
        k = k(ii);
        res(ii) = .5*h2 + .5*(1-k.*h2.*tanh(k*h1))./k.*tanh(k*h2);
    end
%     res = res./(.5*(h2.*sech(k2.*h2).^2+tanh(k2.*h2)./k2));
end
