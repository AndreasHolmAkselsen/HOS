%     1 2   3   4  5  6   7  8
h1  = [.5,.5,.53,.5,.5,.5,.6,.6];
h2 = h1-.42
kh1 = [1.8,2.1,2.5,2.5,3.2,4.1,2.5,3];

k1 = kh1./h1;
lam1 = 2*pi./k1

kh2 = [.54,.6,.65,.68,.78,.91,.97,1.1];
k2 = kh2./h2;
lam2 = 2*pi./k2


 
% Numbers from Li's second paper.
clear
L = 15;

% Case A, C, D
f0 = .6:.05:.95;
hd = .55;
hs__hd = 0.36;
k0hd = [1.03,1.15,1.27,1.40,1.55,1.71,1.88,2.06];
k0shs = [0.57,0.62,0.67,0.73,.79,.85,.91,.97];
kd = k0hd./hd
hs = hd.*hs__hd
ks = k0shs./hs;
lam_d = 2*pi./kd 
lam_s = 2*pi./ks

L./lam_d
L./lam_s

3*(hd-hs)

% case B
hd = .75;
hs__hd = 0.53;
f0 = [.55,.7,.8,1.05,1.4];
k0hd = [1.13,1.60,2.00,3.34,5.92];
k0shs = [.76,1.02,1.22,1.86,3.17];

kd = k0hd./hd
hs = hd.*hs__hd
ks = k0shs./hs;
lam_d = 2*pi./kd 
lam_s = 2*pi./ks

L./lam_d
L./lam_s



