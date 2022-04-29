function [t,Y] = RK4(func,tt,Y0)
% tt = [t0,dt,t_end]
dt = tt(2);
t = (tt(1):dt:tt(3))';
nt = length(t);
Y = zeros(length(Y0),nt);
Y(:,1) = Y0;
for i=1:nt-1 
    y_1 = func(t(i),Y(:,i));
    y_2 = func(t(i)+0.5*dt,Y(:,i)+0.5*dt*y_1);
    y_3 = func(t(i)+0.5*dt,(Y(:,i)+0.5*dt*y_2));
    y_4 = func(t(i)+dt,(Y(:,i)+y_3*dt));

    Y(:,i+1) = Y(:,i) + (1/6)*(y_1+2*y_2+2*y_3+y_4)*dt;
    if any(isnan(Y(:,i+1))), break; end
end
Y =  Y.';
end