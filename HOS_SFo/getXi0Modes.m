function Axi0=getXi0Modes(xi,eta,k,M,h,P)

alpha=1/P;
xi2=xi;
if P>1
    for iP=1:P-1
        eta1=eta*(1-alpha*(iP-1));
        eta2=eta*(1-alpha*iP);
        xi2=H2Operator(xi2,eta1,eta2,k,M,h);
    end
else
    eta2=eta;
end
Axi0=getModeAmplitudes(H2Operator(xi2,eta2,zeros(size(eta)),k,M,h));

end

