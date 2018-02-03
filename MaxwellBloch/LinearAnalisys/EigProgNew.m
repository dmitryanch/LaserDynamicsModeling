function X=EigProgNew(dk,r)  %%%функция с учетом сдвига фазы
global sigma gamma a delta q 
k=sqrt(delta/a);
om=-(2*sqrt(delta*a)*dk+a*dk^2)/(1+sigma); %вспомогательный символ для краткости
rm=1-(om^2+2*sqrt(delta*a)*om*dk+a*om*dk^2)/sigma;
Est=sqrt(r-rm);
Pst=(1+1i*(om+2*sqrt(delta*a)*dk+a*dk^2)/(1+sigma))*Est;
for i=1:length(q)
    
        %для 1Д-неустойчивости
%     M=[ -1i*om-1i*a*q(i)^2-2*1i*a*(k+dk)*q(i)-2*1i*a*k*dk-1i*a*dk^2-sigma 0 sigma 0 0;
%         0 1i*om+1i*a*q(i)^2-2*1i*a*(k+dk)*q(i)+2*1i*a*k*dk+1i*a*dk^2-sigma 0 sigma 0;
%         rm 0 -1i*om-1 0 Est;
%         0 rm 0 1i*om-1 Est;
%         -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*Est -gamma/2*Est -gamma];
    
    %%%% для 2Д=неустойчивости
        M=[ -1i*om-1i*a*q(i)^2-2*1i*a*k*dk-1i*a*dk^2-sigma 0 sigma 0 0;
            0 1i*om+1i*a*q(i)^2+2*1i*a*k*dk+1i*a*dk^2-sigma 0 sigma 0;
            rm 0 -1i*om-1 0 Est;
            0 rm 0 1i*om-1 conj(Est);
            -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];

    X(i)=max(real(eig(M)));
end
end