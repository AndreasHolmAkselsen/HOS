function [zz,df,z] = fmap_SCnum(xx,yy,H,theta,xx_b)
    assert(isvector(xx)&&isvector(yy));
    xx = xx(:)'; yy = yy(:);
    assert(yy(2)>yy(1))
    yy = flipud(yy);
    zz = xx + 1i*yy;
    
    % compute df for each corner
    lambda = exp(zz-xx_b);
    df_i = ((lambda+1)./(lambda+(H(2:end)./H(1:end-1)).^(pi/theta))).^(theta/pi); % df ~ 1/tau
    
    xi_cut = 50; % xi_value at which to follow asymptote in step.
    iPlus = real(zz-xx_b) > xi_cut;
    iMinus = real(zz-xx_b) < -xi_cut;
    df_i(iPlus) = 1;
    temp = H(1:end-1)./H(2:end)+0*zz;
    df_i(iMinus) = temp(iMinus);
   
    df = prod(df_i,3);
    [ny,nx] = size(zz);
    df_xh = .5*(df(1,1:nx-1)+df(1,2:nx));
    df_yh = .5*(df(1:ny-1,:)+df(2:ny,:));
    z = cumsum([0,df_xh.*diff(xx)],2) + cumsum([zeros(1,nx);df_yh.*diff(1i*yy)],1);
%     z0 = interp2(xx,yy,real(z),0,0)+1i*interp2(xx,yy,imag(z),0,0);
    z0 = interp2(xx,yy,real(z),0,-pi)+1i*interp2(xx,yy,imag(z),0,0);
    z = z-z0; % orientation constant
    
    if isscalar(xx_b)
        K = H(2)./pi;
    else
        K = 1/(-min(imag(z(ny,:))))*max(H);% scaling constant
    end
    z = K*z;
    df = K*df;
        
    [zz,df,z] = deal(flipud(zz),flipud(df),flipud(z));
end

% Obs! df = prod(df_i) -> ddf = D[df] ~= prod(ddf_i)!!!
% function [zz,ddf] = ddfz(xx,yy)
%     global H theta xx_b
%     assert(isvector(xx)&&isvector(yy));
%     xx = xx(:)'; yy = yy(:);
%     assert(yy(2)>yy(1))
%     yy = flipud(yy);
%     zz = xx + 1i*yy;
%     
%     % compute df for each corner
%     c = (H(2:end)./H(1:end-1)).^(pi/theta/2);
%     lambda = exp(zz-xx_b);
%     tau_i = ((lambda+c.^2)./(lambda+1)).^(theta/pi); % df ~ 1/tau
%     
%     ddf_i = tau_i.*lambda.*(c.^2-1)./(c.^2+lambda).^2;
%     xi_cut = 50; % xi_value at which to follow asymptote in step.
%     ddf_i( real(zz-xx_b) > xi_cut |  real(zz-xx_b) < -xi_cut ) = 0;
%     ddf = prod(ddf_i,3);
%     [zz,ddf] = deal(flipud(zz),flipud(ddf));    
% end
