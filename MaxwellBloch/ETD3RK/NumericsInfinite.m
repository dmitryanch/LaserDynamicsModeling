clear 
clc

a=sym('a', 'positive');
sigma=sym('sigma', 'positive');
gamma=sym('gamma', 'positive');
r=sym('r', 'positive');
delta=sym('delta', 'real');
K=sym('K', 'positive');
tau=sym('tau','positive');
Ne=3;

M=[-sigma-1i*a*K,sigma,0;
            0,-1-1i*delta,0;
            0,0,-gamma];
        
J1=M;
pairs=[[1,1];[1,2];[2,2];[3,3]];

title='function [';
expM=simplify(expm(J1*tau),'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('expM%d%d = expM(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('expM(pairs(i,1),pairs(i,2))'));
    disp([sprintf('expM%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',expM%d%d', pairs(i,1),pairs(i,2))];
end
exp2M=simplify(expm(J1*tau/2),'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('exp2M%d%d = exp2M(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('exp2M(pairs(i,1),pairs(i,2))'));
    disp([sprintf('exp2M%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',exp2M%d%d', pairs(i,1),pairs(i,2))];
end
invM=simplify(J1^-1,'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('invM%d%d = invM(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('invM(pairs(i,1),pairs(i,2))'));
    disp([sprintf('invM%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',invM%d%d', pairs(i,1),pairs(i,2))];
end

%%%

a2=simplify((exp2M-eye(Ne))*invM,'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('a2%d%d = a2(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('a2(pairs(i,1),pairs(i,2))'));
    disp([sprintf('a2%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',a2%d%d', pairs(i,1),pairs(i,2))];
end
b2=simplify((expM-eye(Ne))*invM,'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('b2%d%d = b2(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('b2(pairs(i,1),pairs(i,2))'));
    disp([sprintf('b2%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',b2%d%d', pairs(i,1),pairs(i,2))];
end
c1=simplify(-4/tau^2*invM^3 - 1/tau*invM^2 + expM*4/tau^2*invM^3 - expM*3/tau*invM^2 + expM*invM,'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('c1%d%d = c1(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('c1(pairs(i,1),pairs(i,2))'));
    disp([sprintf('c1%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',c1%d%d', pairs(i,1),pairs(i,2))];
end
c2=simplify(2/tau^2*invM^3 + 1/tau*invM^2 + expM*(-2)/tau^2*invM^3 + expM/tau*invM^2,'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('c2%d%d = c2(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('c2(pairs(i,1),pairs(i,2))'));
    disp([sprintf('c2%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',c2%d%d', pairs(i,1),pairs(i,2))];
end
c3=simplify(-4/tau^2*invM^3 - 3/tau*invM^2 - invM + expM*4/tau^2*invM^3 - expM/tau*invM^2,'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('c3%d%d = c3(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('c3(pairs(i,1),pairs(i,2))'));
    disp([sprintf('c3%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',c3%d%d', pairs(i,1),pairs(i,2))];
end
disp([title,']=precompBuf()']);
%%%
% expM11=exp(- sigma*tau - K*a*tau*1i);
% expM12=-(sigma*(exp(- sigma*tau - K*a*tau*1i) - exp(- tau - delta*tau*1i))*1i)/(delta + sigma*1i - K*a - 1i);
% expM22=exp(-tau*(delta*1i + 1));
% expM33=exp(-gamma*tau);
% exp2M11=exp(-(sigma*tau)/2)*exp(-(K*a*tau*1i)/2);
% exp2M12=-(sigma*(exp(- (sigma*tau)/2 - (K*a*tau*1i)/2) - exp(- tau/2 - (delta*tau*1i)/2))*1i)/(delta + sigma*1i - K*a - 1i);
% exp2M22=exp(- tau/2 - (delta*tau*1i)/2);
% exp2M33=exp(-(gamma*tau)/2);
% invM11=-(sigma - K*a*1i)/(K^2*a^2 + sigma^2);
% invM12=-sigma/((delta*1i + 1)*(sigma + K*a*1i));
% invM22=(delta*1i - 1)/(delta^2 + 1);
% invM33=-1/gamma;
% a211=-((exp(-(tau*(sigma + K*a*1i))/2) - 1)*(sigma - K*a*1i))/(K^2*a^2 + sigma^2);
% a212=- (sigma*(exp(-(tau*(sigma + K*a*1i))/2) - 1))/((delta*1i + 1)*(sigma + K*a*1i)) - (sigma*(delta*1i - 1)*(exp(-(tau*(sigma + K*a*1i))/2) - exp(-(tau*(delta*1i + 1))/2))*1i)/((delta^2 + 1)*(delta + sigma*1i - K*a - 1i));
% a222=((delta*1i - 1)*(exp(-(tau*(delta*1i + 1))/2) - 1))/(delta^2 + 1);
% a233=-(exp(-(gamma*tau)/2) - 1)/gamma;
% b211=-((exp(-tau*(sigma + K*a*1i)) - 1)*(sigma - K*a*1i))/(K^2*a^2 + sigma^2);
% b212=- (sigma*(exp(-tau*(sigma + K*a*1i)) - 1))/((delta*1i + 1)*(sigma + K*a*1i)) - (sigma*(delta*1i - 1)*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*1i)/((delta^2 + 1)*(delta + sigma*1i - K*a - 1i));
% b222=((delta*1i - 1)*(exp(-tau*(delta*1i + 1)) - 1))/(delta^2 + 1);
% b233=-(exp(-gamma*tau) - 1)/gamma;
% c111=-(exp(-tau*(sigma + K*a*1i))*(exp(tau*(sigma + K*a*1i))*4i - sigma^2*tau^2*1i - sigma*tau*3i - sigma*tau*exp(tau*(sigma + K*a*1i))*1i + 3*K*a*tau + K^2*a^2*tau^2*1i + 2*K*a*sigma*tau^2 + K*a*tau*exp(tau*(sigma + K*a*1i)) - 4i)*1i)/(tau^2*(sigma + K*a*1i)^3);
% c112=(sigma*(delta - sigma*1i + K*a - 1i)*1i)/(tau*(delta - 1i)^2*(sigma + K*a*1i)^2) - (sigma*exp(-tau*(sigma + K*a*1i)))/((delta*1i + 1)*(sigma + K*a*1i)) - (sigma*(delta*1i - 1)*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*1i)/((delta^2 + 1)*(delta + sigma*1i - K*a - 1i)) + (sigma*exp(-tau*(sigma + K*a*1i))*(delta - sigma*1i + K*a - 1i)*3i)/(tau*(delta - 1i)^2*(sigma + K*a*1i)^2) + (sigma*(delta*1i - 1)^2*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*3i)/(tau*(delta^2 + 1)^2*(delta + sigma*1i - K*a - 1i)) - (sigma*(delta*1i - 1)^3*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*4i)/(tau^2*(delta^2 + 1)^3*(delta + sigma*1i - K*a - 1i)) + (8*sigma*sinh((tau*(sigma + K*a*1i))/2)*exp(-(tau*(sigma + K*a*1i))/2)*(delta*2i + sigma + K*a*1i + delta*sigma*1i - delta^2 - K^2*a^2 + sigma^2 - K*a*delta + K*a*sigma*2i + 1))/(tau^2*(sigma*1i - K*a)^3*(delta - 1i)^3);
% c122=(exp(-tau*(delta*1i + 1))*(delta*1i - 1))/(delta^2 + 1) - (4*(delta*1i - 1)^3)/(tau^2*(delta^2 + 1)^3) - (delta*1i - 1)^2/(tau*(delta^2 + 1)^2) - (3*exp(-tau*(delta*1i + 1))*(delta*1i - 1)^2)/(tau*(delta^2 + 1)^2) + (4*exp(-tau*(delta*1i + 1))*(delta*1i - 1)^3)/(tau^2*(delta^2 + 1)^3);
% c133=- exp(-gamma*tau)/gamma - (4*exp(-gamma*tau) + gamma*tau*(3*exp(-gamma*tau) + 1) - 4)/(gamma^3*tau^2);
% c211=(sigma - K*a*1i)^2/(tau*(K^2*a^2 + sigma^2)^2) - (2*(sigma - K*a*1i)^3)/(tau^2*(K^2*a^2 + sigma^2)^3) + (exp(- sigma*tau - K*a*tau*1i)*(sigma - K*a*1i)^2)/(tau*(K^2*a^2 + sigma^2)^2) + (2*exp(- sigma*tau - K*a*tau*1i)*(sigma - K*a*1i)^3)/(tau^2*(K^2*a^2 + sigma^2)^3);
% c212=(sigma*(delta*1i - 1)^3*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*2i)/(tau^2*(delta^2 + 1)^3*(delta + sigma*1i - K*a - 1i)) - (sigma*(delta*1i - 1)^2*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*1i)/(tau*(delta^2 + 1)^2*(delta + sigma*1i - K*a - 1i)) - (4*sigma*sinh((tau*(sigma + K*a*1i))/2)*exp(-(tau*(sigma + K*a*1i))/2)*(delta*2i + sigma + K*a*1i + delta*sigma*1i - delta^2 - K^2*a^2 + sigma^2 - K*a*delta + K*a*sigma*2i + 1))/(tau^2*(sigma*1i - K*a)^3*(delta - 1i)^3) - (sigma*cosh((tau*(sigma + K*a*1i))/2)*exp(-(tau*(sigma + K*a*1i))/2)*(delta - sigma*1i + K*a - 1i)*2i)/(tau*(delta - 1i)^2*(sigma + K*a*1i)^2);
% c222=(delta*1i - 1)^2/(tau*(delta^2 + 1)^2) + (2*(delta*1i - 1)^3)/(tau^2*(delta^2 + 1)^3) + (exp(-tau*(delta*1i + 1))*(delta*1i - 1)^2)/(tau*(delta^2 + 1)^2) - (2*exp(-tau*(delta*1i + 1))*(delta*1i - 1)^3)/(tau^2*(delta^2 + 1)^3);
% c233=(2*exp(-gamma*tau) + gamma*tau*(exp(-gamma*tau) + 1) - 2)/(gamma^3*tau^2);
% c311=(sigma - K*a*1i)/(K^2*a^2 + sigma^2) - (3*(sigma - K*a*1i)^2)/(tau*(K^2*a^2 + sigma^2)^2) + (4*(sigma - K*a*1i)^3)/(tau^2*(K^2*a^2 + sigma^2)^3) - (exp(-tau*(sigma + K*a*1i))*(sigma - K*a*1i)^2)/(tau*(K^2*a^2 + sigma^2)^2) - (4*exp(-tau*(sigma + K*a*1i))*(sigma - K*a*1i)^3)/(tau^2*(K^2*a^2 + sigma^2)^3);
% c312=sigma/((delta*1i + 1)*(sigma + K*a*1i)) + (sigma*(delta - sigma*1i + K*a - 1i)*3i)/(tau*(delta - 1i)^2*(sigma + K*a*1i)^2) + (sigma*exp(-tau*(sigma + K*a*1i))*(delta - sigma*1i + K*a - 1i)*1i)/(tau*(delta - 1i)^2*(sigma + K*a*1i)^2) + (sigma*(delta*1i - 1)^2*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*1i)/(tau*(delta^2 + 1)^2*(delta + sigma*1i - K*a - 1i)) - (sigma*(delta*1i - 1)^3*(exp(-tau*(sigma + K*a*1i)) - exp(-tau*(delta*1i + 1)))*4i)/(tau^2*(delta^2 + 1)^3*(delta + sigma*1i - K*a - 1i)) + (8*sigma*sinh((tau*(sigma + K*a*1i))/2)*exp(-(tau*(sigma + K*a*1i))/2)*(delta*2i + sigma + K*a*1i + delta*sigma*1i - delta^2 - K^2*a^2 + sigma^2 - K*a*delta + K*a*sigma*2i + 1))/(tau^2*(sigma*1i - K*a)^3*(delta - 1i)^3);
% c322=(4*exp(-tau*(delta*1i + 1))*(delta*1i - 1)^3)/(tau^2*(delta^2 + 1)^3) - (3*(delta*1i - 1)^2)/(tau*(delta^2 + 1)^2) - (4*(delta*1i - 1)^3)/(tau^2*(delta^2 + 1)^3) - (exp(-tau*(delta*1i + 1))*(delta*1i - 1)^2)/(tau*(delta^2 + 1)^2) - (delta*1i - 1)/(delta^2 + 1);
% c333=1/gamma - (4*exp(-gamma*tau) + gamma*tau*(exp(-gamma*tau) + 3) - 4)/(gamma^3*tau^2);
% function [,expM11,expM12,expM22,expM33,exp2M11,exp2M12,exp2M22,exp2M33,invM11,invM12,invM22,invM33,a211,a212,a222,a233,b211,b212,b222,b233,c111,c112,c122,c133,c211,c212,c222,c233,c311,c312,c322,c333]=precompBuf()