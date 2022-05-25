function a = icosfft(A)
% Reconstrict a function a(x) from cosine decomposition A(k) along first dimension; 
% a = icosfft(cosfft(a)).
% Equivalent to 
% a = sum([.5*A(1);A(2:end-1);.5*A(end)].'.*cos(k.*x),2)./(n-1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n    = size(A,1);
aPer = ifft([A;flipud(A(2:n-1,:))],[],1);
a    = aPer(1:n,:);

% % test: 
% assert(iscolumn(A));
% a = sum([.5*A(1);A(2:n-1);.5*A(end)].'.*cos((0:n-1).*(0:n-1).'*pi/(n-1)),2)./(n-1);
end




% % alternative function definition where a = sum(A.'.*cos(k.*x),2)
% function a = icosfft(A)
% % equivalent to
% % a = sum(A.'.*cos(k.*x),2)
% % with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.
% 
% n    = size(A,1);
% A([1,n],:) = 2*A([1,n],:);
% aPer = ifft([A;flipud(A(2:n-1,:))],[],1)*(n-1);
% a    = aPer(1:n,:);
% 
% % % test: 
% % assert(iscolumn(A));
% % a = sum(A.'.*cos((0:n-1).*(0:n-1).'*pi/(n-1)),2);
% end