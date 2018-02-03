function [ef,pf,df,ee]=nonlinear1dnondiag(e,p,d,sigma,delta,gamma,r)

global delta_aperture_ON sigma_aperture_ON

    ee=ifft(e);
%     e1=ee(Nr/2,Nr/2);
    pp=ifft(p);
    dd=real(ifft(d));

%     if sigma_aperture_ON==1
%         ef=fft(sigma.*(pp-ee)+ee);
%     else
%         ef=fft(sigma*pp);
%     end
    if delta_aperture_ON==1
        pf=fft(dd.*ee-1i*delta.*pp);
    else
        pf=fft(dd.*ee);
    end
    df=fft(-gamma*(-r+real(ee).*real(pp)+imag(ee).*imag(pp)));
    df=makeReal1d(df);
    
    
%     ef=filter1d(e);
    ef = zeros(size(e));
    pf=filter1d(pf);
    df=filter1d(df);
    
end