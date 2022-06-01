clear

N_wavelenghtsInL = 20;
T = 1;


Tramp = N_wavelenghtsInL*T;
% dt = .1;
% t=0:dt:Tramp;
f0 = 0.01;
shapeFactor = -2*atanh(2*f0-1); % f at t=0
% f = .5*(1+tanh( shapeFactor*(t/Tramp-.5)));
df0 = .5*shapeFactor*sech(.5*shapeFactor)^2/Tramp;
t_add = f0/df0;
% t_foot = -t_add:dt:0;
% f_foot = f0+df0*t_foot;

dt = .1;
t=0:dt:2*Tramp;
f_ramp = @(t) (t>t_add).*.5.*(1+tanh(shapeFactor*((t-t_add)/Tramp-.5))) + df0.*t.*(t<=t_add);
df_ramp = @(t) (t>t_add).*.5.*shapeFactor.*sech(shapeFactor*((t-t_add)/Tramp-.5)).^2/Tramp + df0.*(t<=t_add);
figure,
subplot(211);plot(t, f_ramp(t),[t_add,t_add+Tramp],[.5,.5],'-+k')
subplot(212);plot(t, df_ramp(t),[t_add,t_add+Tramp],[0,0],'-+k')
