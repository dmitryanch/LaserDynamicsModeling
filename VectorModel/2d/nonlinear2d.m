function [x, ep, em, ffttime]=nonlinear2d(X)
global Nr
tic
    ep = ifft2(reshape(X(1,:,:),Nr,Nr));
    em = ifft2(reshape(X(2,:,:),Nr,Nr));
    pp=ifft2(reshape(X(3,:,:),Nr,Nr));
    pm=ifft2(reshape(X(4,:,:),Nr,Nr));
    np=ifft2(reshape(X(5,:,:),Nr,Nr));
    nm=ifft2(reshape(X(6,:,:),Nr,Nr));
    m=ifft2(reshape(X(7,:,:),Nr,Nr));
ffttime = toc;
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
    
tic
    ppf=fft2(pp1);
    pmf=fft2(pm1);
    npf=fft2(np1);
    nmf=fft2(nm1);
    mf=fft2(m1);
ffttime = ffttime + toc;
    
    npf=makeReal2d(npf);
    nmf=makeReal2d(nmf);
    
%     epf=filter2d(epf);
%     emf=filter2d(emf);
    ppf=filter2d(ppf);
    pmf=filter2d(pmf);
    npf=filter2d(npf);
    nmf=filter2d(nmf);
    mf=filter2d(mf);
    
    x=zeros(size(X));
    x(1,:,:)=epf;
    x(2,:,:)=emf;
    x(3,:,:)=ppf;
    x(4,:,:)=pmf;
    x(5,:,:)=npf;
    x(6,:,:)=nmf;
    x(7,:,:)=mf;
    
end