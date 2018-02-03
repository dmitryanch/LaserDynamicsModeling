function [x, ep, em]=nonlinear1d(X)

    ep = ifft(X(1,:));
    em = ifft(X(2,:));
    pp=ifft(X(3,:));
    pm=ifft(X(4,:));
    np=ifft(X(5,:));
    nm=ifft(X(6,:));
    m=ifft(X(7,:));
    
    epf=zeros(size(X(1,:)));
    emf=zeros(size(X(2,:)));
    ppf=fft(-np.*ep-m.*em);
    pmf=fft(-nm.*em-conj(m).*ep);
    npf=fft(1/2*(conj(ep).*pp+ep.*conj(pp))+1/4*(conj(em).*pm+em.*conj(pm)));
    nmf=fft(1/2*(conj(em).*pm+em.*conj(pm))+1/4*(conj(ep).*pp+ep.*conj(pp)));
    mf=fft(1/4*(ep.*conj(pm)+conj(em).*pp));
    
    npf=makeReal1d(npf);
    nmf=makeReal1d(nmf);
    
%     epf=filter1d(epf);
%     emf=filter1d(emf);
    ppf=filter1d(ppf);
    pmf=filter1d(pmf);
    npf=filter1d(npf);
    nmf=filter1d(nmf);
    mf=filter1d(mf);
    
    x=[epf;emf;ppf;pmf;npf;nmf;mf];
end