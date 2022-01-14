function f_x = HOSODE45(t,y) 
global HOSODEsyst
[phiS_t,eta_t] = HOSODEsyst(t,y(1:end/2),y(end/2+1:end));
f_x = [phiS_t;eta_t];
end