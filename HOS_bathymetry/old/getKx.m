function kx = getKx(x,nk)
% nx = size(x,1);
nx = length(x);
if nargin == 1
    nk = nx;
end
dk = 2*pi/(nx*(x(2)-x(1)));
if mod(nk,2)==0
    kx = [0:nk/2-1,-nk/2:-1]'*dk;
else
    kx = [0:(nk-1)/2, -(nk-1)/2:-1]'*dk;
end

end