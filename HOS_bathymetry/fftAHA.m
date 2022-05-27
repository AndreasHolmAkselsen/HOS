function fftAHA(a)

if nargin==0
    rng(2)
    a = rand(11,1);
    % a = tanh(2*x/L)-.4*cos(6*x*L*pi).*exp(x/L)+1;
end


A = cosfft_cleansum(a);
a_ = icosfft_cleansum(A);

% a([1,end],:)=0; % odd function!
% A = sinfft(a);
% a_ = isinfft(A);

figure('color','w');
n = length(a);
plot(0:n-1,a,0:n-1,a_,'--','linewidth',1.5)

end


function A = cosfft(a)
% equivalent to 
% A = sum([a(1);2*a(2:end-1);a(end)].*cos(k.*x),1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n = size(a,1);
A = fft([a;flipud(a(2:n-1,:))],[],1);
A = real(A(1:n,:));

% % test: 
% assert(iscolumn(a));
% A = sum([a(1);2*a(2:end-1);a(end)].*cos((0:n-1).*(0:n-1).'*pi/(n-1)),1).';
end

function a = icosfft(A)
% equivalent to
% a = sum([.5*A(1);A(2:end-1);.5*A(end)].'.*cos(k.*x),2)./(n-1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n    = size(A,1);
aPer = ifft([A;flipud(A(2:n-1,:))],[],1);
a    = aPer(1:n,:);

% % test: 
% assert(iscolumn(A));
% a = sum([.5*A(1);A(2:n-1);.5*A(end)].'.*cos((0:n-1).*(0:n-1).'*pi/(n-1)),2)./(n-1);
end

function A=sinfft(a)
% equivalent to 
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

function a=isinfft(A)
% equivalent to
% sum([.5*A(1);A(2:n-1);.5*A(end)].'.*sin(k.*x),2)./(n-1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n    = size(A,1);
aPer = ifft([A;-flipud(A(2:n-1,:))],[],1);
a    = imag(aPer(1:n,:));

% % test: 
% assert(iscolumn(A));
% a = sum([.5*A(1);A(2:n-1);.5*A(end)].'.*sin((0:n-1).*(0:n-1).'*pi/(n-1)),2)./(n-1);
end






function A = cosfft_cleansum(a)
% equivalent to 
% A = sum([a(1);2*a(2:end-1);a(end)].*cos(k.*x),1)./(n-1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n = size(a,1);
A = fft([a;flipud(a(2:n-1,:))],[],1)/(n-1);
A = real(A(1:n,:));
A([1,n],:) = .5*A([1,n],:);

% test: 
% assert(iscolumn(a));
% A = sum([a(1);2*a(2:end-1);a(end)].*cos((0:n-1).*(0:n-1).'*pi/(n-1)),1).'./(n-1);
% A([1,n],:) = .5*A([1,n],:);
end

function a = icosfft_cleansum(A)
% equivalent to
% a = sum(A.'.*cos(k.*x),2)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n    = size(A,1);
A([1,n],:) = 2*A([1,n],:);
aPer = ifft([A;flipud(A(2:n-1,:))],[],1)*(n-1); % = .5*fft([A;flipud(A(2:n-1,:))],[],1);
a    = aPer(1:n,:);

% % test: 
% assert(iscolumn(A));
% a = sum(A.'.*cos((0:n-1).*(0:n-1).'*pi/(n-1)),2);
end

function A=sinfft_cleansum(a)
% equivalent to 
% A = sum([a(1);2*a(2:end-1);a(end)].*sin(k.*x),1).'/(n-1)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

assert(all([a(1,:);a(end,:)]==0,'all'),'a(1,:) = a(end,:) = 0 assumed since sine is an odd function.')
n = size(a,1);
A = fft([a;-flipud(a(2:n-1,:))],[],1)/(n-1);
A = -imag(A(1:n,:));
A([1,n],:) = .5*A([1,n],:);

% % test: 
% assert(iscolumn(a));
% A = sum([a(1);2*a(2:end-1);a(end)].*sin((0:n-1).*(0:n-1).'*pi/(n-1)),1).'/(n-1);
% A([1,n],:) = .5*A([1,n],:);
end

function a=isinfft_cleansum(A)
% equivalent to
% sum(A.'.*sin(k.*x),2)
% with dx=L/(n-1), x=(0:n-1).'*dx, k=(0:n-1)*pi/L.

n    = size(A,1);
aPer = ifft([A;-flipud(A(2:n-1,:))],[],1)*(n-1);
a    = imag(aPer(1:n,:));

% % test: 
% assert(iscolumn(A));
% a = sum(A.'.*sin((0:n-1).*(0:n-1).'*pi/(n-1)),2);
end