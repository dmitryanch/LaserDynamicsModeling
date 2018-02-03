clear
clc

sigma=sym('sigma','positive');
a=sym('a','positive');
k=sym('k','real');
r=sym('r','positive');
delta=sym('delta','positive');
gamma=sym('gamma','positive');

M=[sigma+1i*a*k^2,-sigma,0;
-r,1+1i*delta,0;
0,0,gamma];

e=eig(-M);

%%
exp1=(a^2*k^4 - 2*a*delta*k^2 + delta^2 + sigma^2 + 2*sigma + 1)/(sigma^2 + 2*sigma + 1);
exp2=-(- a^2*k^4*sigma^2 - a^2*k^4 + 2*a*delta*k^2*sigma^2 + 2*a*delta*k^2 - delta^2*sigma^2 - delta^2 + sigma^4 + 2*sigma^3 + 2*sigma^2 + 2*sigma + 1)/(sigma*(2*sigma^2 + 4*sigma + 2));

sigma0=0.1;
a0=0.01;
delta0=-3;
gamma0=0.01;
k0=sqrt(delta/a);
K0=-20:0.1:20;

plot(K0,subs(subs(subs(subs(subs(exp1,sigma,sigma0),delta,delta0),a,a0),gamma,gamma0),k,K0),'red');             %% first eigen value
hold on;
plot(K0,subs(subs(subs(subs(subs(exp2,sigma,sigma0),delta,delta0),a,a0),gamma,gamma0),k,K0),'green');           %% second eigen value

plot(K0,1+((a0*K0.^2-delta0)./(1+sigma0)).^2,'blue');                                                           %% expression from D-analisys. It is the same as first eigen value
