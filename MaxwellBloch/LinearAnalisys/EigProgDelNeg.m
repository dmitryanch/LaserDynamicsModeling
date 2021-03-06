% function X=EigProgDelNeg(delta,r)  %%%������� ��� ��������� �����������
function [X Y]=EigProgDelNeg(delta,r)
global sigma gamma a q 

rm=1+delta^2/(1+sigma)^2;
Est=sqrt(r-rm);
Pst=Est*(1-1i*delta/(1+sigma));
for i=1:length(q)
    
        M=[ -(1i*a*q(i)^2+sigma-1i*sigma*delta/(1+sigma)) 0 sigma 0 0;
            0 1i*a*q(i)^2-sigma-1i*sigma*delta/(1+sigma) 0 sigma 0;
            rm 0 -1-1i*delta+1i*sigma*delta/(1+sigma) 0 Est;
            0 rm 0 -1+1i*delta-1i*sigma*delta/(1+sigma) conj(Est);
            -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];

        S=eig(M);
    X(i)=max(real(S));
%     X(:,i)=sort(real(eig(M)));
%     X(i)=max(eig(M));
    
    for i1=1:5
        if real(S(i1))==X(i);
            Y(i)=S(i1);
            break;
        end;
    end
end


%��������� ���� ����������� ����� ��� ���� q
% figure;[X Y]=EigProgDelNeg(delta,r);plot(q,X);%title(['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)]);X(:,1)
%���������� ���� ���������� ����������� ������������ �����
% i1=0;for i=2:(numel(q)-1),if X(5,i-1)<X(5,i) && X(5,i)>X(5,i+1),i1=i1+1;disp(['qmax',num2str(i1),'=',num2str(q(i)),' i=',num2str(i),' Increment=',num2str(X(5,i))]);end,end
% i1=0;for i=2:(numel(q)-1),if X(i-1)<X(i) && X(i)>X(i+1),i1=i1+1;disp(['qmax',num2str(i1),'=',num2str(q(i)),' i=',num2str(i),' Increment=',num2str(Y(i))]);end,end
