function [X Y]=EigProgDelNegRandK(delta,r)  %%%функция для статейных результатов
global sigma gamma a q k
% k=0;
exp1=(a*k^2-delta)/(1+sigma);
exp2=(sigma*delta+a*k^2)/(1+sigma);
rm=1+exp1^2;
Est=sqrt(r-rm);
Pst=Est*(1+1i*(a*k^2-delta)/(1+sigma));
for i=1:length(q)
    
%         M=[ -(1i*a*q(i)^2+2*1i*a*q(i)*k+1i*a*k^2+sigma-1i*exp2) 0 sigma 0 0;
%             0 1i*a*q(i)^2-2*1i*a*q(i)*k+1i*a*k^2-sigma-1i*exp2 0 sigma 0;
%             rm 0 -1-1i*delta+1i*exp2 0 Est;
%             0 rm 0 -1+1i*delta-1i*exp2 conj(Est);
%             -gamma/2*conj(Pst) -gamma/2*Pst -gamma/2*conj(Est) -gamma/2*Est -gamma];

        M=[-sigma -exp2+a*(q(i)+k)^2 0 0 sigma 0 0 0 0 0;
           exp2-a*(q(i)+k)^2 -sigma 0 0 0 sigma 0 0 0 0;
           0 0 -sigma -exp2+a*(q(i)-k)^2 0 0 sigma 0  0 0;
           0 0 exp2-a*(q(i)-k)^2 -sigma 0 0 0 sigma 0 0;
           rm 0 0 0 -1 -exp2+delta 0 0 Est 0;
           0 rm 0 0 exp2-delta -1 0 0 0 Est;
           0 0 rm 0 0 0 -1 -exp2+delta Est 0;
           0 0 0 rm 0 0 exp2-delta -1 0 Est;
           -gamma/2*Est -gamma/2*Est*exp1 -gamma/2*Est -gamma/2*Est*exp1 -gamma/2*Est 0 -gamma/2*Est 0 -gamma 0;
           gamma/2*Est*exp1 -gamma/2*Est -gamma/2*Est*exp1 gamma/2*Est 0 -gamma/2*Est 0 gamma/2*Est 0 -gamma];
        
    S=eig(M);
        [X(i) ind]=max(real(S));
%         Y(i)=S(ind);
        realS=real(S(ind));
        ss=sort(real(S));
        if abs(ss(1)-ss(2))<1e-10
            imagS=abs(imag(S(ind)));
        else 
            imagS=imag(S(ind));
        end
        Y(i)=realS+1i*imagS;
%     X(:,i)=sort(real(S));
%      Y(:,i)=S;
%     for i1=1:numel(S)
%         if real(S(i1))==X(i);
%             Y(i)=S(i1);
%             break;
%         end;
%     end
end
% Y=1;
end

%Диаграмма всех собственных чисел для всех q
% figure;X=EigProgDelNegRandK(delta,r);plot(q,X);title(['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)]);X(:,1)
%Нахождение всех максимумов наибольшего собственного числа
% i1=0;for i=2:(numel(q)-1),if X(5,i-1)<X(5,i) && X(5,i)>X(5,i+1),i1=i1+1;disp(['qmax',num2str(i1),'=',num2str(q(i)),' Increment=',num2str(X(5,i))]);end,end
