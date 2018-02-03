function x=makeReal2d(f)
    global Nr kfil
    
    f(1,1)=real(f(1,1));
%     ! 2 - 4 quarters without axis values
    for i=2:kfil
      for j=2:kfil
        f(Nr+2-i,Nr+2-j)=conj(f(i,j));
      end
    end
%     ! 1 - 3 quarters without axis values
    for i=2:kfil
      for j=Nr-kfil+2:Nr
        f(Nr+2-i,Nr+2-j)=conj(f(i,j));
      end
    end
%     ! x-axis
      for i=2:kfil
        f(Nr+2-i,1)=conj(f(i,1));
        f(1,Nr+2-i)=conj(f(1,i));
      end

x=f;
end