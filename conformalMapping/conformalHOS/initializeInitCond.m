function [eta,H] = initializeInitCond(x_target,eta_target,h_target,nIt)
xi = x_target;
eta = eta_target;
H = h_target;
kx = getKx(xi);


% figure('color','w');hold on
% clf; hold on
% plot(x_target,eta_target,'k.','DisplayName','Target')

x = x_target;
nx = numel(x);
L = (x(2)-x(1))*nx; assert(all(diff(x,2)==0));
% H = min(H,realmax);
for i = 1:nIt
    FFTeta = fft(eta);
    if isfinite(H)
        Lsin = -2./(exp(2*kx.*H)-1-2*(kx==0));
        f = xi + 1i*ifft(FFTeta.*Lsin,[],1);
    else
        f =  xi + 2i*fft(conj(FFTeta/nx).*(kx>0),[],1)+1i*FFTeta(1)/nx;
    end
        
    x_new = real(f);
%     plot(x_new,eta,'.','DisplayName',"Iteration "+i)
    
    eta = interp1([x-L;x;x+L],[eta;eta;eta],x_new,'linear',nan);
    x = x_new;
    H = h_target+mean(eta);
end
% legend

end
