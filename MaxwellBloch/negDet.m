clear  
clc

a=sym('a', 'positive');
delta=sym('delta', 'real');
gamma=sym('gamma', 'positive');
sigma=sym('sigma', 'positive');
r=sym('r', 'positive');
q=sym('q', 'positive');
k=sym('k', 'real');
x=sym('x', 'positive');
t=sym('t', 'real');
lyambda=sym('lyambda', 'real');

syms e1 e2 p1 p2 d

omega=delta*sigma/(sigma+1)
Dst=1+(delta/(sigma+1))^2 
Est=sqrt(r-Dst)%*exp(-1i*omega*t)
Pst=(sigma-1i*omega)*Est/sigma
%%
% D=Dst + d*exp(1i*q*x+lyambda*t)+conj(d)*exp(-1i*q*x+conj(lyambda)*t)
% E=(abs(Est)+e1*exp(1i*q*x+lyambda*t)+e2*exp(-1i*q*x+conj(lyambda)*t))*exp(-1i*omega*t)
% P=(abs(Pst)+p1*exp(1i*q*x+lyambda*t)+p2*exp(-1i*q*x+conj(lyambda)*t))*exp(-1i*omega*t)
% 
% f1=-1i*a*k^2+sigma*(P-E);
% f2=-(1+1i*delta)*P+D*E;
% f3=-gamma*(D-r+0.5*(conj(E)*P+E*conj(P)));
% 
% f1=simplify(f1)
% f2=simplify(f2)
% f3=simplify(f3)

M=[ -(1i*a*q^2+sigma-1i*sigma*delta/(1+sigma)) 0 sigma 0 0;
            0 1i*a*q^2-sigma-1i*sigma*delta/(1+sigma) 0 sigma 0;
            Dst 0 -1-1i*delta+1i*sigma*delta/(1+sigma) 0 Est;
            0 Dst 0 -1+1i*delta-1i*sigma*delta/(1+sigma) conj(Est);
            -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];