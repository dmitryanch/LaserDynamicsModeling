function X=EigProgDelPol(delta,r)  %%%������� ��� ��������� �����������
global sigma gamma a q 

rm=1;
Est=sqrt(r-rm);
Pst=Est;
for i=1:length(q)
    
        %��� 1�-��������������
%     M=[ -(1i*(2*sqrt(a*delta)*q(i)+a*q(i)^2)+sigma) 0 sigma 0 0;
%         0 -(1i*(2*sqrt(a*delta)*q(i)-a*q(i)^2)+sigma) 0 sigma 0;
%         rm 0 -1 0 Est;
%         0 rm 0 -1 conj(Est);
%         -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];
    
    %%%% ��� 2�=��������������
        M=[ -(1i*a*q(i)^2+sigma) 0 sigma 0 0;
            0 -(-1i*a*q(i)^2+sigma) 0 sigma 0;
            rm 0 -1 0 Est;
            0 rm 0 -1 conj(Est);
            -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];

    X(i)=max(real(eig(M)));
end
end