clear all
clc

a=sym('a', 'positive');
b=sym('b', 'positive');
c=sym('c', 'positive');
omega=sym('omega', 'real');
K=sym('K', 'positive');
tau=sym('tau','positive');
gammaA=sym('gammaA','positive');
gammaP=sym('gammaP','positive');
Ne=7;

M=[-1i*a*K-1,-(gammaA+1i*gammaP),0,0,0,0,0;
            -(gammaA+1i*gammaP),-1i*a*K-1,0,0,0,0,0;
            0,0,-1-1i*omega,0,0,0,0;
            0,0,0,-1-1i*omega,0,0,0;
            0,0,0,0,-b,0,0;
            0,0,0,0,0,-b,0;
            0,0,0,0,0,0,-c];
        
J1=M;
pairs=[[1,1];[1,2];[2,1];[2,2];[3,3];[4,4];[5,5];[6,6];[7,7]];

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

a2=simplify((exp2M-eye(7))*invM,'Steps',200);
for i=1:size(pairs,1)
%     eval(sprintf('a2%d%d = a2(pairs(i,1),pairs(i,2));', pairs(i,1),pairs(i,2)));
    expr=char(eval('a2(pairs(i,1),pairs(i,2))'));
    disp([sprintf('a2%d%d=', pairs(i,1),pairs(i,2)),expr,';']);
    title=[title,sprintf(',a2%d%d', pairs(i,1),pairs(i,2))];
end
b2=simplify((expM-eye(7))*invM,'Steps',200);
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
% expM11=cosh(gammaA*tau + gammaP*tau*1i)*exp(- tau - K*a*tau*1i);
% expM12=-sinh(gammaA*tau + gammaP*tau*1i)*exp(- tau - K*a*tau*1i);
% expM21=-sinh(gammaA*tau + gammaP*tau*1i)*exp(- tau - K*a*tau*1i);
% expM22=cosh(gammaA*tau + gammaP*tau*1i)*exp(- tau - K*a*tau*1i);
% expM33=exp(-tau*(omega*1i + 1));
% expM44=exp(-tau*(omega*1i + 1));
% expM55=exp(-b*tau);
% expM66=exp(-b*tau);
% expM77=exp(-c*tau);
% exp2M11=cosh((gammaA*tau)/2 + (gammaP*tau*1i)/2)*exp(- tau/2 - (K*a*tau*1i)/2);
% exp2M12=-sinh((gammaA*tau)/2 + (gammaP*tau*1i)/2)*exp(- tau/2 - (K*a*tau*1i)/2);
% exp2M21=-sinh((gammaA*tau)/2 + (gammaP*tau*1i)/2)*exp(- tau/2 - (K*a*tau*1i)/2);
% exp2M22=cosh((gammaA*tau)/2 + (gammaP*tau*1i)/2)*exp(- tau/2 - (K*a*tau*1i)/2);
% exp2M33=exp(- tau/2 - (omega*tau*1i)/2);
% exp2M44=exp(- tau/2 - (omega*tau*1i)/2);
% exp2M55=exp(-(b*tau)/2);
% exp2M66=exp(-(b*tau)/2);
% exp2M77=exp(-(c*tau)/2);
% invM11=-(K*a*1i + 1)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% invM12=(gammaA + gammaP*1i)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% invM21=(gammaA + gammaP*1i)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% invM22=-(K*a*1i + 1)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% invM33=(omega*1i - 1)/(omega^2 + 1);
% invM44=(omega*1i - 1)/(omega^2 + 1);
% invM55=-1/b;
% invM66=-1/b;
% invM77=-1/c;
% a211=- ((exp(-(tau*(K*a*1i + 1))/2)*cosh((tau*(gammaA + gammaP*1i))/2) - 1)*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (exp(-(tau*(K*a*1i + 1))/2)*sinh((tau*(gammaA + gammaP*1i))/2)*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% a212=((exp(-(tau*(K*a*1i + 1))/2)*cosh((tau*(gammaA + gammaP*1i))/2) - 1)*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (exp(-(tau*(K*a*1i + 1))/2)*sinh((tau*(gammaA + gammaP*1i))/2)*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% a221=((exp(-(tau*(K*a*1i + 1))/2)*cosh((tau*(gammaA + gammaP*1i))/2) - 1)*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (exp(-(tau*(K*a*1i + 1))/2)*sinh((tau*(gammaA + gammaP*1i))/2)*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% a222=- ((exp(-(tau*(K*a*1i + 1))/2)*cosh((tau*(gammaA + gammaP*1i))/2) - 1)*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (exp(-(tau*(K*a*1i + 1))/2)*sinh((tau*(gammaA + gammaP*1i))/2)*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% a233=((omega*1i - 1)*(exp(-(tau*(omega*1i + 1))/2) - 1))/(omega^2 + 1);
% a244=((omega*1i - 1)*(exp(-(tau*(omega*1i + 1))/2) - 1))/(omega^2 + 1);
% a255=-(exp(-(b*tau)/2) - 1)/b;
% a266=-(exp(-(b*tau)/2) - 1)/b;
% a277=-(exp(-(c*tau)/2) - 1)/c;
% b211=- ((exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i)) - 1)*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% b212=((exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i)) - 1)*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% b221=((exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i)) - 1)*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% b222=- ((exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i)) - 1)*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1);
% b233=((omega*1i - 1)*(exp(-tau*(omega*1i + 1)) - 1))/(omega^2 + 1);
% b244=((omega*1i - 1)*(exp(-tau*(omega*1i + 1)) - 1))/(omega^2 + 1);
% b255=-(exp(-b*tau) - 1)/b;
% b266=-(exp(-b*tau) - 1)/b;
% b277=-(exp(-c*tau) - 1)/c;
% c111=(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1)/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (12*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (3*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (12*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (4*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (6*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c112=(exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (4*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (2*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2) + (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (3*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (12*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (4*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (6*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c121=(exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (4*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (2*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2) + (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (3*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (12*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (4*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (6*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c122=(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1)/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (12*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1))/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (3*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (12*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (4*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (6*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c133=(exp(-tau*(omega*1i + 1))*(omega*1i - 1))/(omega^2 + 1) - (4*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) - (omega*1i - 1)^2/(tau*(omega^2 + 1)^2) - (3*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) + (4*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3);
% c144=(exp(-tau*(omega*1i + 1))*(omega*1i - 1))/(omega^2 + 1) - (4*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) - (omega*1i - 1)^2/(tau*(omega^2 + 1)^2) - (3*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) + (4*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3);
% c155=- exp(-b*tau)/b - (4*exp(-b*tau) + b*tau*(3*exp(-b*tau) + 1) - 4)/(b^3*tau^2);
% c166=- exp(-b*tau)/b - (4*exp(-b*tau) + b*tau*(3*exp(-b*tau) + 1) - 4)/(b^3*tau^2);
% c177=- exp(-c*tau)/c - (4*exp(-c*tau) + c*tau*(3*exp(-c*tau) + 1) - 4)/(c^3*tau^2);
% c211=(6*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (6*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1)/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (2*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (2*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c212=(2*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2) + (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (6*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c221=(2*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2) + (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (6*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c222=(6*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (6*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1)/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (2*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (2*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c233=(omega*1i - 1)^2/(tau*(omega^2 + 1)^2) + (2*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) + (exp(-tau*(omega*1i + 1))*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) - (2*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3);
% c244=(omega*1i - 1)^2/(tau*(omega^2 + 1)^2) + (2*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) + (exp(-tau*(omega*1i + 1))*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) - (2*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3);
% c255=(2*exp(-b*tau) + b*tau*(exp(-b*tau) + 1) - 2)/(b^3*tau^2);
% c266=(2*exp(-b*tau) + b*tau*(exp(-b*tau) + 1) - 2)/(b^3*tau^2);
% c277=(2*exp(-c*tau) + c*tau*(exp(-c*tau) + 1) - 2)/(c^3*tau^2);
% c311=(K*a*1i + 1)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (3*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (12*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (12*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (4*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c312=(6*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2) - (4*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (gammaA + gammaP*1i)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (12*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (4*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (2*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c321=(6*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2) - (4*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (gammaA + gammaP*1i)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) - (exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (12*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (4*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (2*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c322=(K*a*1i + 1)/(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1) + (3*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) + (12*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) + (exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*2i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - K^2*a^2 + 1))/(tau*(2*K*a - 2*gammaA*gammaP + gammaA^2*1i - gammaP^2*1i + K^2*a^2*1i - 1i)^2) - (12*exp(-tau*(K*a*1i + 1))*cosh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*((K*a*2i)/3 + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - (K^2*a^2)/3 + 1/3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (4*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(gammaA + gammaP*1i)*(K*a*6i + gammaA*gammaP*2i + gammaA^2 - gammaP^2 - 3*K^2*a^2 + 3))/(tau^2*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^3) - (2*exp(-tau*(K*a*1i + 1))*sinh(tau*(gammaA + gammaP*1i))*(K*a*1i + 1)*(gammaA + gammaP*1i))/(tau*(K*a*2i - gammaA*gammaP*2i - gammaA^2 + gammaP^2 - K^2*a^2 + 1)^2);
% c333=(4*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) - (3*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) - (4*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) - (exp(-tau*(omega*1i + 1))*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) - (omega*1i - 1)/(omega^2 + 1);
% c344=(4*exp(-tau*(omega*1i + 1))*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) - (3*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) - (4*(omega*1i - 1)^3)/(tau^2*(omega^2 + 1)^3) - (exp(-tau*(omega*1i + 1))*(omega*1i - 1)^2)/(tau*(omega^2 + 1)^2) - (omega*1i - 1)/(omega^2 + 1);
% c355=1/b - (4*exp(-b*tau) + b*tau*(exp(-b*tau) + 3) - 4)/(b^3*tau^2);
% c366=1/b - (4*exp(-b*tau) + b*tau*(exp(-b*tau) + 3) - 4)/(b^3*tau^2);
% c377=1/c - (4*exp(-c*tau) + c*tau*(exp(-c*tau) + 3) - 4)/(c^3*tau^2);
% function [,expM11,expM12,expM21,expM22,expM33,expM44,expM55,expM66,expM77,exp2M11,exp2M12,exp2M21,exp2M22,exp2M33,exp2M44,exp2M55,exp2M66,exp2M77,invM11,invM12,invM21,invM22,invM33,invM44,invM55,invM66,invM77,a211,a212,a221,a222,a233,a244,a255,a266,a277,b211,b212,b221,b222,b233,b244,b255,b266,b277,c111,c112,c121,c122,c133,c144,c155,c166,c177,c211,c212,c221,c222,c233,c244,c255,c266,c277,c311,c312,c321,c322,c333,c344,c355,c366,c377]=precompBuf()