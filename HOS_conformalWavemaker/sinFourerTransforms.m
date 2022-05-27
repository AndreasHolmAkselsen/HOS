clear
% clc
 
rng(2)

L = 1;
Nx = 11;
dx = L/(Nx-1);
x = (0:Nx-1).'*dx;
k = (0:Nx-1)*pi/L; 

a = rand(Nx,1);
% a = tanh(2*x/L)-.4*cos(6*x*L*pi).*exp(x/L)+1;
a([1,end])=0; % odd function!

x2 = (0:2*Nx-1)'*dx; % to plot more of the domain
% x2 = x;

% a=a-mean(a);


%% using seb's. Fully fft
% a2=a.';
% aPer=[a2,-fliplr(a2(2:end-1))];
% A=fft(aPer,[],2);
% A=-imag(A(:,1:Nx));
% aPer=ifft([A,-fliplr(A(2:end-1))],[],2);
% a_=imag(aPer(:,1:Nx));

%% transfrom fft
% a2=a.';
% aPer=[a2,-fliplr(a2(2:end-1))];
% A=fft(aPer,[],2);
% A=-imag(A(:,1:Nx));
% a_ = sum([A(1),2*A(2:end-1),A(end)].*sin(k.*x2),2)./(2*(Nx-1)) ;

%% transfrom fft, pure sum definition|||
% a2=a.';
% aPer=[a2,-fliplr(a2(2:end-1))];
% A=fft(aPer,[],2);
% A=-imag(A(:,1:Nx))./(2*(Nx-1)); A(2:end-1) = 2*A(2:end-1);
% a_ = sum(A.*sin(k.*x2),2) ;

%% inverse transfrom fft
% A = sum([a(1);2*a(2:end-1);a(end)].*sin(k.*x),1);
% aPer=ifft([A,-fliplr(A(2:end-1))],[],2);
% a_=imag(aPer(:,1:Nx));

%% fully sum
A = sum([a(1);2*a(2:end-1);a(end)].*sin(k.*x),1);
a_ = sum([.5*A(1),A(2:end-1),.5*A(end)].*sin(k.*x2),2)./(Nx-1) ;
 


fprintf('%.3f,\t',a);fprintf('\n');
fprintf('%.3f,\t',a_);fprintf('\n');

plot(x,a,x2,a_,'--')
