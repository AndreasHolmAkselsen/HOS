function A=sinfft(a)
% Sine decomposition of function a(x) along first dimension; a = isinfft(sinfft(a)).
% Equivalent to 
% A = sum([a(1);2*a(2:end-1);a(end)].*sin(k.*x),1).'
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

assert(all([a(1,:);a(end,:)]==0,'all'),'a(1,:) = a(end,:) = 0 assumed since sine is an odd function.')
n = size(a,1);
A = fft([a;-flipud(a(2:n-1,:))],[],1);
A = -imag(A(1:n,:));

% % test: 
% assert(iscolumn(a));
% A = sum([a(1);2*a(2:end-1);a(end)].*sin((0:n-1).*(0:n-1).'*pi/(n-1)),1).';
end



% % alternative function definition where a = sum(A.'.*sin(k.*x),2)
% function A=sinfft(a)
% % equivalent to 
% % A = sum([a(1);2*a(2:end-1);a(end)].*sin(k.*x),1).'/(n-1)
% % with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.
% 
% assert(all([a(1,:);a(end,:)]==0,'all'),'a(1,:) = a(end,:) = 0 assumed since sine is an odd function.')
% n = size(a,1);
% A = fft([a;-flipud(a(2:n-1,:))],[],1)/(n-1);
% A = -imag(A(1:n,:));
% A([1,n],:) = .5*A([1,n],:);
% 
% % % test: 
% % assert(iscolumn(a));
% % A = sum([a(1);2*a(2:end-1);a(end)].*sin((0:n-1).*(0:n-1).'*pi/(n-1)),1).'/(n-1);
% % A([1,n],:) = .5*A([1,n],:);
% end