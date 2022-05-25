function k = findWaveNumbers(omega,h,Ux,ni)
% k = findWaveNumbers(omega,h,Ux,ni) vecotrized solver for the generalized dispersion relation 
%       (omega - k*Ux)^2 = k*g*tanh(k*h) for real and imaginary solutions of k.
%
%   INPUT: 
%         omega:   [1,nj,...]   - angular frequency [rad/s]
%         h:       scalar       - water depth [m]
%         Ux:      scalar       - current in wave propactaion direction (=U cos(theta) where theta is the angle between current direction and wave vector).
%         ni:      scalar       - number of imaginary roots to be found. (ni==0 gives real root only.)
%
%   OUTPUT:
%         k        [ni,nj,...]  - roots of the dispersion relation.
%                                 k(1,:) contians the real root of the respective frequencies of sign omega.
%                                 k(2,:end,:) contains the negative imaginary roots (1i*k_i>0) of the respective frequencies.
%                                 If Ux==0, l = -k will be the complimentary roots of the dispersion relation.
%                                 (Note that there may exist up to four real roots for nonzero Ux!)
%
%   AHA, 03/07-2020

    if numel(omega)==1 && omega==0 && ni==0, k=0; return; end

%     if size(omega,1)==1, omega = repmat(omega,[ni+1,1]); end
    
    
    g = 9.81;
    c = omega*sqrt(h/g);
    if nargin<3 || Ux == 0
        x = xTanhx_minus_cSq_Newton(c,ni);
    else
        d = Ux/sqrt(g*h);
        x = xTanhx_minus_cxdSq_Newton(c,d,ni);
    end
    k = x/h;
    
    k(1,omega==0) = 0;
%     k(repmat(omega==0,[ni+1,1]))==0
end


function x = xTanhx_minus_cSq_Newton(c,ni)
% solves x.*tanh(x) - c.^2 = 0 for x finding one real and ni imaginary solutions.
    csq = c.^2;
%     ii = [-ni:-1,1:ni]'; % x_{-i} will equal -x_i
%     ii = (1:ni)'; % the positive imaginary solutions
    ii = -(1:ni)';
    x = [c.*abs(c)  %csq % for sign(omega) use c.*abs(c)  
        1i*(ii*pi-sign(ii).*atan(csq))]; % initial guess (first element for the real solution, the rest for the imaginaty)
    
%     x = [c(1,:,:,:,:,:).*abs(c(1,:,:,:,:,:))  %csq % for sign(omega) use c.*abs(c)  
%         1i*(ii*pi-sign(ii).*atan(csq(1,:,:,:,:,:)))]; % initial guess (first element for the real solution, the rest for the imaginaty)
    
    
    N = 100;
    error_threshold = 1e-9;
    for n = 1:N
       x = x - ( x.*tanh(x)-csq )./( x.*sech(x).^2+tanh(x) ); 
       err = abs((x.*tanh(x)-csq)); % or  std(...
       max_err = max(err(:));
       if max_err <  error_threshold
           break
       end
    end
    if n == N, fprintf('failed to converge withing n=%d, largest errro = %g.\n',N,max_err); end
    
%   %   validataion plot    
%     hf=figure;hf.Position(3)=800;
%     xx = linspace(0,5,1000);
%     subplot(1,2,1);
%     plot(xx,xx.*tanh(xx)-csq,'k-', xx,xx-csq,'k--', x(1),0,'ko',xx,0*xx,'k:')
%     ylim([-10,10])
%     subplot(1,2,2);
%     xx = linspace(ii(1)*pi-pi,ii(end)*pi,1000);
%     yy = 1i*xx;
%     plot(xx,yy.*tanh(yy)-csq,'k-', xx,csq+sign(xx).*tan(xx),'k--', imag(x(2:end)),0,'ko',xx,0*xx,'k:')
%     ylim([-10,10])
end



function x = xTanhx_minus_cxdSq_Newton(c,d,ni)
% solves x.*tanh(x) - (c-d.*x).^2 = 0 for x finding one real and ni imaginary solutions.
    ii = -(1:ni)';
    x = [abs(.5*(1+2*c.*d-sqrt(1+4*c.*d))/d.^2).*sign(c) % for sign omega, use abs(.5*(1+2*c.*d-sqrt(1+4*c.*d))/d.^2).*sign(c)
        1i*(ii*pi-sign(ii).*atan(c.^2))];  % initial guess (first element for the real solution, the rest for the imaginaty)    
    N = 100;
    error_threshold = 1e-9;
    for n = 1:N
        x = x - ( x.*tanh(x)-(c-d.*x).^2 )./( 2*d.*(x-d.*x) + x.*sech(x).^2+tanh(x) );
        err = abs((x.*tanh(x)-(c-d.*x).^2 )); % or  std(...
        max_err = max(err(:));
        if max_err <  error_threshold
            break
        end
    end
    if n == N, fprintf('failed to converge withing n=%d, largest errro = %g.\n',N,max_err); end
end