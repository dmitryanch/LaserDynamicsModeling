function X=EigProg(k,r)  %%%функция для статейных результатов
global sigma gamma a delta q 
rm=1+(a*k^2-delta)^2/(1+sigma)^2;
z=(sigma*delta+a*k^2)/(1+sigma); %вспомогательный символ для краткости
Est=sqrt(r-rm);
Pst=(1+1i*(a*k^2-delta)/(1+sigma))*Est;
for i=1:length(q)
    
        %для 1Д-неустойчивости
    M=[ -sigma+1i*(z-a*(k+q(i))^2) 0 sigma 0 0;
        0 -sigma-1i*(z-a*(k-q(i))^2) 0 sigma 0;
        rm 0 1i*(z-delta)-1 0 Est;
        0 rm 0 -1i*(z-delta)-1 conj(Est);
        -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];
    
    %%%% для 2Д=неустойчивости
%         M=[-sigma+1i*(z-a*(k^2+q(i)^2)) 0 sigma 0 0;
%             0 -sigma-1i*(z-a*(k^2+q(i)^2)) 0 sigma 0;
%             rm 0 1i*(z-delta)-1 0 Est;
%             0 rm 0 -1i*(z-delta)-1 conj(Est);
%             -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];

    X(i)=max(real(eig(M)));
%     X(:,i)=sort(real(eig(M)));
end
end