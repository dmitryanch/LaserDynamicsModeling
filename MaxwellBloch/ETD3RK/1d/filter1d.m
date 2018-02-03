function x=filter1d(f)
global Nr kfil
    x=f;
    x(kfil+1:Nr-kfil+1)=0;
end