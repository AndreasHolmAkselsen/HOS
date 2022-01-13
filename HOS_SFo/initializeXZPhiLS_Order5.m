function [deltax,deltaz,phiLS]=initializeXZPhiLS_Order5(k,a,h,alpha)

%Miscellaneous coefficients
A13=1/2;
C13=(1/2-A13)/3; %From mass conservation

%c-value at free surface 
c=0;  

%Dispersion relation 
g=9.81;
T=tanh(k*h);
sig0=sqrt(g*k*T);

%Horizontal particle grid
Psi=k*alpha;

%First-order solution
x1=-a*sin(Psi)*exp(k*c);
z1=a*cos(Psi)*exp(k*c);
phiLS1=a*sig0/k*sin(Psi)*exp(k*c);

%Second-order solution
x2=zeros(size(Psi));
z2=ones(size(Psi))*0.5*k*a^2*exp(2*k*c);

%Third-order solution
x3=k^2*a^3*(-2*exp(3*k*c)+(1/2)*exp(k*c))*sin(Psi);
z3=k^2*a^3*(exp(3*k*c)-(1/2)*exp(k*c))*cos(Psi);

%Fourth-order solution
x4=k^3*a^4*((1/6)*exp(4*k*c)-(1/2)*exp(2*k*c))*sin(2*Psi);
z4=k^3*a^4*((-(1/3)*exp(4*k*c)+(1/2)*exp(2*k*c))*cos(2*Psi)+((3/2)*exp(4*k*c)-(1/2)*exp(2*k*c)));
phiLS4=a^4*sig0*k^2*((-(2/3)*exp(4*k*c)+(1/2)*exp(2*k*c))*sin(2*Psi));

%Fifth-order solution
x5=k^4*a^5*((-53/8*exp(5*k*c)+A13*exp(3*k*c)+7/6*exp(k*c))*sin(Psi)+(1/72*exp(5*k*c)-1/12*exp(3*k*c))*sin(3*Psi));
z5=k^4*a^5*((21/8*exp(5*k*c)+C13*exp(3*k*c)-7/6*exp(k*c))*cos(Psi)+(-1/24*exp(5*k*c)+1/12*exp(3*k*c))*cos(3*Psi));
phiLS5=a^5*sig0*k^3*((-3/8*exp(5*k*c)+1/2*exp(3*k*c)-19/24*exp(k*c))*sin(Psi)+(-1/8*exp(5*k*c)+1/12*exp(3*k*c))*sin(3*Psi));

%Return particle displacement and surface potential
deltax=(x1+x2+x3+x4+x5);
deltaz=(z1+z2+z3+z4+z5);
phiLS=(phiLS1+phiLS4+phiLS5);

end