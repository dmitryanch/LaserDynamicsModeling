clear
clear v q del Ne
tic
syms v q del Ne
A=[0 q^2 v/2;
   q^2 0 -v*del/2;
   2*(1-Ne/(1+del^2)) 0 (1-Ne/(1+del^2))];
A2=subs(A,[v q],[20 (2*Ne/abs(del))^0.5]); %подстановка определенных ню v и ку максимальное
c=-0.1:-0.025:-4.5; %диапазон изменени€ и шаг по дельта
b=0;
for n=1:size(c,2)
    m=1+c(n)^2;o=1; %лева€ граница диапазона Ne
    while m<=20 %права€ граница Ne
        S=subs(A2,{del Ne},{c(n) m});
        L(:,n,o)=eig(S); %вектор собственных значений
        h=0;
        for f=1:size(L,1)
          if real(L(f,n,o))>0
             h=h+1;
          end;
        end;
        if h>=2 %условие по количеству положительных реальных частей собственных значений
            b=b+1; %счетчик по точкам из области неустойчивости
            delta(b)=c(n);N(b)=m; %точки области неустойчивости
        end;
        m=m+0.05;o=o+1;
    end;
end;

%---------------- ¬ыделение границ области неустойчивости-----------------

p=0;
for i=1:b
    if delta(i)==delta(1)
        p=p+1;delta1(p)=delta(i);N1(p)=N(i);
    elseif delta(i)<delta(i-1) && N(i-1)+0.05<20
        p=p+1;delta1(p)=delta(i-1);N1(p)=N(i-1);
    end;
end;
for i=1:p
    del1(i)=delta1(p+1-i);Ne1(i)=N1(p+1-i);
end;
for i=2:b
    if delta(i)<delta(i-1)
        del1(size(del1,2)+1)=delta(i);Ne1(size(Ne1,2)+1)=N(i);
    end;
end;

%-------------------------------------------------------------------------

toc
plot(Ne1,-del1);
xlabel('Ne');ylabel('-delta');title('v=20');grid on
%figure,plot(N,-delta,1+c.^2,-c);
%xlabel('Ne');ylabel('-delta');title('v=10,√=0.02');grid on