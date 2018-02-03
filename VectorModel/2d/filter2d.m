function x=filter2d(f)
global Nr kfil
    if numel(size(f)) >2 
        x=reshape(f,Nr,Nr);
    else
        x=f;
    end
    for i=1:Nr
        x(i,kfil+1:Nr-kfil+1)=0;
        x(kfil+1:Nr-kfil+1,i)=0;
    end
end