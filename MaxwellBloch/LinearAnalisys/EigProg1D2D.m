function X=EigProg1D2D(k,r)  %%%функция для программы с учетом 4х направлений одновременно на основе статейного подхода
global sigma gamma a delta q 
rm=1+(a*k^2-delta)^2/(1+sigma)^2;
w=(sigma*delta+a*k^2)/(1+sigma); %вспомогательный символ для краткости (-омега)
Est=sqrt(r-rm);
Pst=(1+1i*(a*k^2-delta)/(1+sigma))*Est;
for i=1:length(q)
    
        M=[ 1i*(w-a*(q(i)+k)^2)-sigma 0 0 0 sigma 0 0 0 0 0;
            0 -1i*(w-a*(q(i)-k)^2)-sigma 0 0 0 sigma 0 0 0 0;
            0 0 1i*(w-a*(q(i)^2+k^2))-sigma 0 0 0 sigma 0 0 0;
            0 0 0 -1i*(w-a*(q(i)^2+k^2))-sigma 0 0 0 sigma 0 0;
            rm 0 0 0 1i*w-1-1i*delta 0 0 0 Est 0;
            0 rm 0 0 0 -1i*w-1+1i*delta 0 0 Est 0;
            0 0 rm 0 0 0 1i*w-1-1i*delta 0 0 Est;
            0 0 0 rm 0 0 0 -1i*w-1+1i*delta 0 Est;
            -gamma/2*conj(Pst) -gamma/2*Pst 0 0 -gamma/2*Est -gamma/2*Est 0 0 -gamma 0;
            0 0 -gamma/2*conj(Pst) -gamma/2*Pst 0 0 -gamma/2*Est -gamma/2*Est 0 -gamma];

    X(i)=max(real(eig(M)));
end
end