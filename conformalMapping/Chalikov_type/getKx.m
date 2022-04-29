function kx = getKx(x)

% nx = size(x,1);
nx = length(x);
dk = 2*pi/(nx*(x(2)-x(1)));
if mod(nx,2)==0
    kx = [0:nx/2-1,-nx/2:-1]'*dk;
else
    kx = [0:(nx-1)/2, -(nx-1)/2:-1]'*dk;
end

end