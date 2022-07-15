function [f,X] = FFTfreq(dt,x,orientation)
% no orientation input gives f=[0,+df,...,-N*df,...,-df]
% orientation = 'centre' gives f=[-...,-df,0,+df,...]
% orientation = 'positive' gives only positive reqencies f=[0,+df,...]
n = size(x,1);
X = fft(x,[],1);
df=1/(n*dt);

if mod(n,2)==0
    f = [0:n/2-1,-n/2:-1]'*df;
else
    f = [0:(n-1)/2, -(n-1)/2:-1]'*df;   
end

if nargin>2 
    switch orientation
        case 'centre'
            f = fftshift(f); X = fftshift(X);
        case 'positive'
            X(f<0,:) = [];
            f(f<0) = [];
        otherwise
            error('Orientation input "%s" not recongnised.',orientation)
    end
end

end