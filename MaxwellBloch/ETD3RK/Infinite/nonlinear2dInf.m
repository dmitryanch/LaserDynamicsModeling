function [ef,pf,df,ee]=nonlinear2dInf(e,p,d,delta,gamma,r)

global delta_aperture_ON Nr

    ee=ifft2(e);
%     e1=ee(Nr/2,Nr/2);
    pp=ifft2(p);
    dd=real(ifft2(d));

    if delta_aperture_ON==1
        pf=fft2(dd.*ee-1i*delta.*pp);
    else
        pf=fft2(dd.*ee);
    end
    df=fft2(-gamma*(-r+real(ee).*real(pp)+imag(ee).*imag(pp)));
    df=makeReal2d(df);
    
    
    ef=zeros(Nr,Nr);
    pf=filter2d(pf);
    df=filter2d(df);
    
end