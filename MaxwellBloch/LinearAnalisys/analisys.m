% clc 
clear
global sigma gamma a q k
q=0:0.01:100;
k=0;
gamma=5e-4;
sigma=3e-2;
delta=-0.5;    
a=0.01;
r=2;
% k=0;
% gamma=0.2;
% sigma=0.01;
% delta=-3;    
% a=1e-4;
% r=25;
Dst=1+((a*k^2+delta)/(1+sigma));
Ist=r-Dst;
Qunst = sqrt(sqrt(2*gamma*sigma*Ist/(1+delta^2)/a^2));
disp(Qunst);
[x y] = EigProgDelNegRandK(delta,r);
% [x1 y1] = EigProgDelNeg(delta,r);
plot(q,x,'-k');
% figure;plot(q,x1);