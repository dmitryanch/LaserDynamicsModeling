function X=EigProgSVL(k,r)  %%%функция для статейных результатов
global sigma gamma a delta q 
rm=1+(a*k^2-delta)^2/(1+sigma)^2;
z=(sigma*delta+a*k^2)/(1+sigma); %вспомогательный символ для краткости
Est=sqrt(r-rm);
Pst=(1+1i*(a*k^2-delta)/(1+sigma))*Est;
for i=1:length(q)
    
        %для 1Д-неустойчивости
    M=[ ];

    E1=-1i*omega*e1-1i*a*(k+q(i))^2*e1+sigma*(p1-e1);
    E2=-1i*omega*e2-1i*a*(k-q(i))^2*e2+sigma*(p2-e2);
    E3=-1i*omega*e3-1i*a*(k-q(i))^2*e3+sigma*(p3-e3);
    E4=-1i*omega*e4-1i*a*(k+q(i))^2*e4+sigma*(p4-e4);
    E5=-1i*omega*e5-1i*a*(k+q(i))^2*e5+sigma*(p5-e5);
    E6=-1i*omega*e6-1i*a*(k-q(i))^2*e6+sigma*(p6-e6);
    E7=-1i*omega*e7-1i*a*(k-q(i))^2*e7+sigma*(p7-e7);
    E8=-1i*omega*e8-1i*a*(k+q(i))^2*e8+sigma*(p8-e8);
    
    P1=-1i*omega*p1-(1+1i*delta)*p1+A*(d1+d3)+e1*B+A^2*e3;
    P2=-1i*omega*p2-(1+1i*delta)*p2+A*(d1'+d4)+e2*B+A^2*e4;
    P3=-1i*omega*p3-(1+1i*delta)*p3+A*(d1+d4')+e3*B+A^2*e1;
    P4=-1i*omega*p4-(1+1i*delta)*p4+A*(d1'+d3')+e4*B+A^2*e2;
    P5=-1i*omega*p5-(1+1i*delta)*p5+A*(d2-d5)+e5*B+A^2*e7;
    P6=-1i*omega*p6-(1+1i*delta)*p6+A*(d2'-d6)+e6*B+A^2*e8;
    P7=-1i*omega*p7-(1+1i*delta)*p7+A*(d6'-d2)+e7*B+A^2*e5;
    P8=-1i*omega*p8-(1+1i*delta)*p8+A*(d5'-d2')+e8*B+A^2*e6;
    
    D1=-gamma(d1+0.5*A*(e3+p2'+e1+p4'+e2'+p1+e4'+p3));
    D2=-gamma(d2+0.5*A*(e5+p6'-e7-p8'+e6'+p5-e8'-p7));
    D3=-gamma(d3+0.5*A*(e1+p4'+p1+e4'));
    D4=-gamma(d4+0.5*A*(e2+p3'+p2+e3'));
    D5=-gamma(d5+0.5*A*(p8'-e5+e8'-p5));
    D6=-gamma(d6+0.5*A*(p7'-e6+e7'-p6));
        
    %%%% для 2Д=неустойчивости
%         M=[-sigma+1i*(z-a*(k^2+q(i)^2)) 0 sigma 0 0;
%             0 -sigma-1i*(z-a*(k^2+q(i)^2)) 0 sigma 0;
%             rm 0 1i*(z-delta)-1 0 Est;
%             0 rm 0 -1i*(z-delta)-1 conj(Est);
%             -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];

    X(i)=max(real(eig(M)));
end
end