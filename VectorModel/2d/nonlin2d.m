function [epf,emf,ppf,pmf,npf,nmf,mf, ep, em]=nonlin2d(Epf,Emf,Ppf,Pmf,Npf,Nmf,Mf)
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
    epf=zeros(size(ep));
    emf=zeros(size(em));
%     ppf=fft2(-np.*ep-m.*em);
%     pmf=fft2(-nm.*em-conj(m).*ep);
%     npf=fft2(1/2*(conj(ep).*pp+ep.*conj(pp))+1/4*(conj(em).*pm+em.*conj(pm)));
%     nmf=fft2(1/2*(conj(em).*pm+em.*conj(pm))+1/4*(conj(ep).*pp+ep.*conj(pp)));
%     mf=fft2(1/4*(ep.*conj(pm)+conj(em).*pp));
    pp1=(-np.*ep-m.*em);
    pm1=(-nm.*em-conj(m).*ep);
    np1=(1/2*(conj(ep).*pp+ep.*conj(pp))+1/4*(conj(em).*pm+em.*conj(pm)));
    nm1=(1/2*(conj(em).*pm+em.*conj(pm))+1/4*(conj(ep).*pp+ep.*conj(pp)));
    m1=(1/4*(ep.*conj(pm)+conj(em).*pp));
    
% tic
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