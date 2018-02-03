function [epf,emf,ppf,pmf,npf,nmf,mf, ep, em]=nonlin1d(Epf,Emf,Ppf,Pmf,Npf,Nmf,Mf,sigma,r)
global Nr
% tic
    ep = ifft(Epf);
    em = ifft(Emf);
    pp=ifft(Ppf);
    pm=ifft(Pmf);
    np=ifft(Npf);
    nm=ifft(Nmf);
    m=ifft(Mf);
% ffttime = toc;
    ep1=sigma.*(pp-ep)+ep;
    em1=sigma.*(pm-em)+em;
    pp1=(-np.*ep-m.*em+r.*ep);
    pm1=(-nm.*em-conj(m).*ep+r.*em);
    np1=(1/2*(conj(ep).*pp+ep.*conj(pp))+1/4*(conj(em).*pm+em.*conj(pm)));
    nm1=(1/2*(conj(em).*pm+em.*conj(pm))+1/4*(conj(ep).*pp+ep.*conj(pp)));
    m1=(1/4*(ep.*conj(pm)+conj(em).*pp));
    
% tic
    epf=fft(ep1);
    emf=fft(em1);
    ppf=fft(pp1);
    pmf=fft(pm1);
    npf=fft(np1);
    nmf=fft(nm1);
    mf=fft(m1);
% ffttime = ffttime + toc;
    
    npf=makeReal1d(npf);
    nmf=makeReal1d(nmf);
    
%     epf=filter1d(epf);
%     emf=filter1d(emf);
    ppf=filter1d(ppf);
    pmf=filter1d(pmf);
    npf=filter1d(npf);
    nmf=filter1d(nmf);
    mf=filter1d(mf);    
end