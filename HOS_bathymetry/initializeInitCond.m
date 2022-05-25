function [eta,xi] = initializeInitCond(fz,x_target,h_target,xi,eta0,nIt)

eta = eta0;


figure('color','w');hold on
clf; hold on
plot(x_target,h_target,'k.','DisplayName','Target')

x = x_target;
nx = numel(x);
L = (x(2)-x(1))*nx; assert(all(abs(diff(x,2))<1e-12*x(end)));

h_ip = griddedInterpolant([x_target-L;x_target;x_target+L],[h_target;h_target;h_target]);
dyy = 1e-4;
for i = 1:nIt      
    
    z = fz(xi+1i*eta);

    dydyy = .5*imag(fz(xi+1i*(eta+dyy))-fz(xi+1i*(eta-dyy)))/dyy;

    
    eta = eta - (h_target-h_ip(real(z)))./dydyy;
    plot(z,'.','DisplayName',"Iteration "+i)
    
end
legend

end
