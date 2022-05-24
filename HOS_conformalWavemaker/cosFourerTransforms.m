clear
% clc
 
rng(2)

L = 1;
Nx = 10;
N = Nx-1;
dx = L/(Nx-1);
x = (0:Nx-1).'*dx;
k = (0:Nx-1)*pi/L; 

a = rand(Nx,1);
% a = tanh(2*x/L)-.4*cos(6*x*L*pi).*exp(x/L)+1;


% x2 = (0:2*Nx-1)'*dx; % to plot more of the domain
x2 = x;


%% using seb's. Fully fft
% x2 = x; % necessary
% a2=a.';
% aPer=[a2,fliplr(a2(:,2:end-1))];
% A=fft(aPer,[],2);
% A=real(A(:,1:size(a2,2)));
% aPer=ifft([A,fliplr(A(:,2:end-1))],[],2);
% a_=aPer(:,1:size(A,2));
 

%% inverse sum transform
% a2=a.';aPer=[a2,fliplr(a2(:,2:end-1))];
% A=fft(aPer,[],2);
% A=real(A(:,1:size(a2,2)));
% % %%
% a_ = sum([A(1),2*A(2:end-1),A(end)].*cos(k.*x2),2)./(2*(Nx-1));


%% fully sum
A = sum([a(1);2*a(2:end-1);a(end)].*cos(k.*x),1);
a_ = sum([A(1),2*A(2:end-1),A(end)].*cos(k.*x2),2)./(2*(Nx-1)) ;


fprintf('%.3f,\t',a);fprintf('\n');
fprintf('%.3f,\t',a_);fprintf('\n');

plot(x,a,x2,a_,'--')
