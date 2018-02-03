function [epf,emf,ppf,pmf,npf,nmf,mf, ep, em]=nonlin2dBuf(Epf,Emf,Ppf,Pmf,Npf,Nmf,Mf,sigma,r)
global Nr
% tic
    ep = ifft2(Epf);
    em = ifft2(Emf);
    pp=ifft2(Ppf);
    pm=ifft2(Pmf);
    np=ifft2(Npf);
    nm=ifft2(Nmf);
    m=ifft2(Mf);
% ffttime = toc;
    ep1=sigma.*(pp-ep)+ep;
    em1=sigma.*(pm-em)+em;
    pp1=(-np.*ep-m.*em+r.*ep);
    pm1=(-nm.*em-conj(m).*ep+r.*em);
    np1=(1/2*(conj(ep).*pp+ep.*conj(pp))+1/4*(conj(em).*pm+em.*conj(pm)));
    nm1=(1/2*(conj(em).*pm+em.*conj(pm))+1/4*(conj(ep).*pp+ep.*conj(pp)));
    m1=(1/4*(ep.*conj(pm)+conj(em).*pp));
    
% tic
    epf=fft2(ep1);
    emf=fft2(em1);
    ppf=fft2(pp1);
    pmf=fft2(pm1);
    npf=fft2(np1);
    nmf=fft2(nm1);
    mf=fft2(m1);
% ffttime = ffttime + toc;
    
    npf=makeReal2d(npf);
    nmf=makeReal2d(nmf);
    
%     epf=filter2d(epf);
%     emf=filter2d(emf);
    ppf=filter2d(ppf);
    pmf=filter2d(pmf);
    npf=filter2d(npf);
    nmf=filter2d(nmf);
    mf=filter2d(mf);    
end