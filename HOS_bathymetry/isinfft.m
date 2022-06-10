function a=isinfft(A)
% Reconstrict a function a(x) from sine decomposition A(k) along first dimension; 
% a = isinfft(sinfft(a)).
% equivalent to
% sum([.5*A(1);A(2:n-1);.5*A(end)].'.*sin(k.*x),2)./(n-1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n    = size(A,1);
aPer = -1i*ifft([A;-flipud(A(2:n-1,:))],[],1);
a    = aPer(1:n,:);

% % test: 
% assert(iscolumn(A));
% a = sum([.5*A(1);A(2:n-1);.5*A(end)].'.*sin((0:n-1).*(0:n-1).'*pi/(n-1)),2)./(n-1);
end




% % alternative function definition where a = sum(A.'.*sin(k.*x),2)
% function a=isinfft(A)
% % equivalent to
% % sum(A.'.*sin(k.*x),2)
% % with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.
% 
% n    = size(A,1);
% aPer = ifft([A;-flipud(A(2:n-1,:))],[],1)*(n-1);
% a    = imag(aPer(1:n,:));
% 
% % % test: 
% % assert(iscolumn(A));
% % a = sum(A.'.*sin((0:n-1).*(0:n-1).'*pi/(n-1)),2);
% end