function [eta,phiS,W,hphi,kx] = HOSODEeqCurr_steadyState(x,df,nIt,relax)
g = 9.81;
[eta,phiS_x,eta_x,W] = deal(0*x);

% figure('color','w','position',[ -1704 505 560 420]); hp=plot(x,0*x,'k'); 
% [WArr, etaArr] = deal(zeros(length(x),nIt));

for it = 1:nIt
    dfS = df(x + 1i*eta);
    Phi_x =  real(dfS);
    Phi_z = -imag(dfS);
    
    
%     WNew =  (phiS_x.*eta_x + Phi_x.*eta_x - Phi_z) - eta_x.^2.*W ;
    WNew =  (phiS_x.*eta_x + Phi_x.*eta_x - Phi_z)./(1+ eta_x.^2)  ;
    etaNew = - ( .5*(phiS_x+Phi_x).^2 +.5*Phi_z.^2  - .5*(1+eta_x.^2).*W.^2)/g; %  - (Phi_x.*eta_x - Phi_z).*W /g ;
    
    maxChangeW = max(abs(WNew-W));
    maxChangeEta = max(abs(etaNew-eta));
    fprintf('iteration %d: maxChangeW = %g, maxChangeEta = %g\n',it, maxChangeW,maxChangeEta);
    
%     W = WNew;
%     eta = etaNew;
    W = relax*WNew + (1-relax)*W;
    eta = relax*etaNew + (1-relax)*eta;
    if maxChangeW < 5e-5 && maxChangeEta < 5e-5, break; end
    
%     WArr(:,it) = W;
%     etaArr(:,it) = eta;
    
    [phiS_x,eta_x] = phiComponentsHOS_steadyState(W,eta);
%     hp.YData = eta; drawnow;
end

% figure();
% subplot(2,1,1); plot(x,WArr); legend("it = "+(1:nIt)');ylabel('W');
% subplot(2,1,2); plot(x,etaArr); legend("it = "+(1:nIt)');ylabel('\eta');
% fprintf('it %d, max change: %.3g\n',[(2:nIt);(max(diff(WArr,1,2),[],1))])
% fprintf('it %d, max change: %.3g\n',[(2:nIt);(max(diff(etaArr,1,2),[],1))])

[~,~,phiS,hphi,kx] = phiComponentsHOS_steadyState(W,eta);

end
