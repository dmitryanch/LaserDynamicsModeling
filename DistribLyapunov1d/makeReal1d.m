function x=makeReal1d(f)
    global Nr kfil
    
    f(1,1)=real(f(1,1));
    for i=2:kfil
      f(Nr+2-i)=conj(f(i));
    end
x=f;
end