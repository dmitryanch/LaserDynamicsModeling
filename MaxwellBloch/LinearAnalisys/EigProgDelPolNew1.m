function X=EigProgDelPolNew1(k,r)  %%%������� ��� ���������� �����: E=Est+e1*exp(iqx+iwt)+e2*exp(-iqx+iwt)
global sigma gamma a q delta
rm=1+(a*k^2-delta)^2/(1+sigma)^2;
z=(sigma*delta+a*k^2)/(1+sigma); %��������������� ������ ��� ���������
Est=sqrt(r-rm);
Pst=(1+1i*(a*k^2-delta)/(1+sigma))*Est;
for i=1:length(q)
    
    M=[ -1i*a*q(i)^2-sigma+1i*z 0 sigma 0 0;
        0 1i*a*q(i)^2-sigma-1i*z 0 sigma 0;
        rm 0 -1-1i*(delta-z) 0 Est;
        0 rm 0 -1+1i*(delta-z) conj(Est);
        -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];

X(i)=max(real(eig(M)));
% X(:,i)=sort(real(eig(M)));

end
end

%��������� ���� ����������� ����� ��� ���� q
% figure;X=EigProgDelNegRandK(delta,r);plot(q,X);title(['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)]);X(:,1)
%���������� ���� ���������� ����������� ������������ �����
% i1=0;for i=2:(numel(q)-1),if X(5,i-1)<X(5,i) && X(5,i)>X(5,i+1),i1=i1+1;disp(['qmax',num2str(i1),'=',num2str(q(i)),' Increment=',num2str(X(5,i))]);end,end
