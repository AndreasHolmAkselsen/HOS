function x = setReal(x,name)
maxRelImag = mean(abs(imag(x)./(abs(x)+1e-9)),'all');
if maxRelImag > 1e-6
    if nargin == 1
        fprintf('Complex values! Mean rel imag = %e\n',maxRelImag);
    else
        fprintf('Complex values in %s! Mean rel imag = %e\n',name,maxRelImag);
    end
end
x = real(x);