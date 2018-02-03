function [Xx KK, Xnorm]=getDistribByWaveNumber2D(absX,kx,N)
global Nr
dk=max(abs(kx))/N;
KK=0:dk:max(abs(kx));
Xx=zeros(size(KK));
Xnorm=zeros(size(KK));
prev=0;
for n=1:numel(KK),    
    limK=KK(n);
    num=0;
    for i=1:Nr
        for j=1:Nr
            if sqrt(kx(i)^2+kx(j)^2)<=limK
                Xx(n)=Xx(n)+absX(i,j);
                num=num+1;
            end
        end
    end
    if n > 1
        Xx(n)=Xx(n)-prev;
    end
    Xnorm(n)=Xx(n)/num;
    prev=prev+Xx(n);
end