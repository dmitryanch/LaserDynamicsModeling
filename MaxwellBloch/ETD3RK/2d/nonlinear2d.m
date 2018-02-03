function [ef,pf,df,ee]=nonlinear2d(e,p,d,sigma,delta,gamma,r)

global delta_aperture_ON sigma_aperture_ON

    ee=ifft2(e);
%     e1=ee(Nr/2,Nr/2);
    pp=ifft2(p);
    dd=real(ifft2(d));

    if sigma_aperture_ON==1
        ef=fft2(sigma.*(pp-ee)+ee);
    else
        ef=fft2(sigma*pp);
    end
    if delta_aperture_ON==1
        pf=fft2(dd.*ee-1i*delta.*pp);
    else
        pf=fft2(dd.*ee);
    end
    df=fft2(-gamma*(-r+real(ee).*real(pp)+imag(ee).*imag(pp)));
    df=makeReal2d(df);
    
    
    ef=filter2d(ef);
    pf=filter2d(pf);
    df=filter2d(df);
    
end