function [ef,pf,df,ee]=nonlinear1dInf(e,p,d,delta,gamma,r)

global delta_aperture_ON Nr

    ee=ifft(e);
%     e1=ee(Nr/2,Nr/2);
    pp=ifft(p);
    dd=real(ifft(d));

    if delta_aperture_ON==1
        pf=fft(dd.*ee-1i*delta.*pp);
    else
        pf=fft(dd.*ee);
    end
    df=fft(-gamma*(-r+real(ee).*real(pp)+imag(ee).*imag(pp)));
    df=makeReal1d(df);
    
    
    ef=zeros(1,Nr);
    pf=filter1d(pf);
    df=filter1d(df);
    
end