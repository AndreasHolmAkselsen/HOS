

fun=@(x) [x(1)^2-exp(x(1)*x(2)),x(1)+x(2).^3];


[x, resnorm, f, exitflag, output, jacob] = newtonraphson(fun, [1,1]);

x
f


[X1,X2] = meshgrid(-1:.05:1);
figure, mesh(X1,X2,X1.^2-exp(X1.*X2));hold on;
mesh(X1,X2,X1+X2.^3);
 plot(x(1),x(2),'*r','markersize',10);
 surf(X1,X2,0*X1)