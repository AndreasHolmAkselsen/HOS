function f_x = HOSODE45(t,y) 
[phiS_t,eta_t] = HOSODEeqCurr(t,y(1:end/2),y(end/2+1:end));
f_x = [phiS_t;eta_t];
end