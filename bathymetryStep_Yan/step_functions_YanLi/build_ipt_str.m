function [ipt] = build_ipt_str(cv,t_vec)
ipt.h_d     = cv(1);
ipt.h_s     = cv(1) -0.35;
ipt.f0      = cv(2);
ipt.omega_0 = 2*pi*cv(2);
ipt.delta   = cv(3);
ipt.phase   = cv(4);
ipt.epsilon = cv(5);
ipt.t_vec   = t_vec;
ipt.x_0     = 0.9;
ipt.T_0     = 32 ;
end

