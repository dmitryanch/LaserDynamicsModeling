function X=EigProgExp(k,r)  %%%функция для статейных результатов
global sigma gamma a delta q 
rm=1+(a*k^2-delta)^2/(1+sigma)^2;
z=(sigma*delta+a*k^2)/(1+sigma); %вспомогательный символ для краткости
Est=sqrt(r-rm);
y=(a*k^2-delta)/(1+sigma);
for i=1:length(q)
    
        %для 1Д-неустойчивости
%     M=[     -sigma -(z-a*(k+q(i))^2) 0 0 sigma 0 0 0 0 0;
%             (z-a*(k+q(i))^2) -sigma 0 0 0 sigma 0 0 0 0;
%             0 0 -sigma -(z-a*(k-q(i))^2) 0 0 sigma 0 0 0;
%             0 0 (z-a*(k-q(i))^2) -sigma 0 0 0 sigma 0 0;
%             rm 0 0 0 -1 delta-z 0 0 Est 0;
%             0 rm 0 0 z-delta -1 0 0 0 Est;
%             0 0 rm 0 0 0 -1 delta-z Est 0;
%             0 0 0 rm 0 0 z-delta -1 0 -Est;
%             -gamma/2*conj(Est) -gamma/2*conj(Est)*y -gamma/2*Est -gamma/2*Est*y -gamma/2*conj(Est) 0 -gamma/2*Est 0 -gamma 0;
%             gamma/2*conj(Est)*y -gamma/2*conj(Est) -gamma/2*Est*y gamma/2*Est 0 -gamma/2*conj(Est) 0 gamma/2*Est 0 -gamma];
    
    %%%% для 2Д=неустойчивости
        M=[     -sigma -(z-a*(k^2+q(i)^2)) 0 0 sigma 0 0 0 0 0;
            (z-a*(k^2+q(i)^2)) -sigma 0 0 0 sigma 0 0 0 0;
            0 0 -sigma -(z-a*(k^2+q(i)^2)) 0 0 sigma 0 0 0;
            0 0 (z-a*(k^2+q(i)^2)) -sigma 0 0 0 sigma 0 0;
            rm 0 0 0 -1 delta-z 0 0 Est 0;
            0 rm 0 0 z-delta -1 0 0 0 Est;
            0 0 rm 0 0 0 -1 delta-z Est 0;
            0 0 0 rm 0 0 z-delta -1 0 -Est;
            -gamma/2*conj(Est) -gamma/2*conj(Est)*y -gamma/2*Est -gamma/2*Est*y -gamma/2*conj(Est) 0 -gamma/2*Est 0 -gamma 0;
            gamma/2*conj(Est)*y -gamma/2*conj(Est) -gamma/2*Est*y gamma/2*Est 0 -gamma/2*conj(Est) 0 gamma/2*Est 0 -gamma];

    X(i)=max(real(eig(M)));
end
end