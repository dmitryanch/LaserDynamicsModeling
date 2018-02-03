function [efr,efi,pfr,pfi,df]=nonlinear1dLyapunov(er,ei,pr,pi,dd,e,p,d,sigma,delta,gamma,r)

global sigma_aperture_ON a

    ep=ifft(e).';
    pp=ifft(p).';
    dp=real(ifft(d)).';

    eer=real(ifft(er));
    eei=real(ifft(ei));
    ppr=real(ifft(pr));
    ppi=real(ifft(pi));
    ddd=real(ifft(dd));

%     if sigma_aperture_ON==1
%         efr=fft(sigma.*(ppr-eer)+eer);
%         efi=fft(sigma.*(ppi-eei)+eei);
%     else
%         efr=fft(sigma*ppr+a*K.*eei);
%         efi=fft(sigma*ppi-a*K.*eer);
%     end
    pfr=fft(dp.*eer+delta.*ppi+real(ep).*ddd);
    pfi=fft(dp.*eei-delta.*ppr+imag(ep).*ddd);
    
    df=fft(-gamma*(real(pp).*eer+imag(pp).*eei+real(ep).*ppr+imag(ep).*ppi));
    
%     efr=filter1d(er);
%     efi=filter1d(ei);
    pfr=filter1d(pfr);
    pfi=filter1d(pfi);
    df=filter1d(df);
    
%     efr=makeReal1d(efr);
%     efi=makeReal1d(efi);
    efr=zeros(size(er));
    efi=zeros(size(ei));
    pfr=makeReal1d(pfr);
    pfi=makeReal1d(pfi);
    df=makeReal1d(df);
end