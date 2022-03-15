function eta = initializeInitCond(x_target,eta_target,H,nIt)
k_cut = inf; % (M+5)*k0;

xi = x_target;
eta = eta_target;
kx = getKx(xi);


% figure('color','w');hold on
% clf; hold on
% plot(x_target,eta_target,'k.','DisplayName','Target')

x = x_target;
nx = numel(x);
L = (x(2)-x(1))*nx; assert(all(diff(x,2)==0));
% argH = exp(-kx.*H)./sinh(kx.*H); argH(1) = 0; argH(kx*H<-10) = -2;
argH = 2./(exp((2*kx.*H))-1); argH(1) = 0;
for i = 1:nIt
    FFTeta = fft(eta);
    if isfinite(H)
        f = xi - 1i*ifft(FFTeta.*argH.*(abs(kx)<k_cut),[],1)+1i*FFTeta(1)/nx;
    else
        f =  xi + 2i*fft(conj(FFTeta/nx).*(abs(kx)<k_cut&kx>0),[],1)+1i*FFTeta(1)/nx;
    end
        
    x_new = real(f);
%     plot(x_new,eta,'.','DisplayName',"Iteration "+i)
    
    eta = interp1([x-L;x;x+L],[eta;eta;eta],x_new);
    x = x_new;
end
% legend

end
