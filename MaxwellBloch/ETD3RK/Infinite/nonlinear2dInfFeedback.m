function [ef,pf,df,ee]=nonlinear2dInfFeedback(e,p,d,delta,gamma,r,sigma,beta)

global delta_aperture_ON Nr

    ee=ifft2(e);
%     e1=ee(Nr/2,Nr/2);
    pp=ifft2(p);
    dd=real(ifft2(d));
    
    ef=fft2(-sigma*beta*abs(ee).^2);
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