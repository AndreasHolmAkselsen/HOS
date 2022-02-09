function [x0,x,z,phiLS,T]=initializeXZPhiLS_RegLagLinear(Lx,Nx,t0,a,k,h)

%Third order dispersion relation at free surface
g=9.81;
T=tanh(k*h);
sig0=sqrt(g*k*T);
om=sig0;
T=2*pi/om;

%Horizontal particle grid
dx=Lx/Nx;
x0=dx*(0:(Nx-1));
Psi=k*x0-om*t0;

%First-order solution
x1=-a*sin(Psi);
z1=a*cos(Psi);
phiLS1=a*sig0/k*sin(Psi);

%Return Lagrangian solution, correct up to order three
x=x0+x1;
z=z1;
phiLS=phiLS1;

end