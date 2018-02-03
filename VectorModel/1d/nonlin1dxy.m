function [exf,eyf,pxf,pyf,nxf,nyf,mf, Ex, Ey]=nonlin1dxy(Exf,Eyf,Pxf,Pyf,Nxf,Nyf,Mf)
global Nr

    Ex = ifft2(Exf);
    Ey = ifft2(Eyf);
    Px=ifft2(Pxf);
    Py=ifft2(Pyf);
    Nx=ifft2(Nxf);
    Ny=ifft2(Nyf);
    M=ifft2(Mf);

px1 = - (Ex*M)/2 - (Ey*M*1i)/2 - (2^(1/2)*Ex*Nx)/2 + (2^(1/2)*Ey*Ny)/2 - (conj(M)*(Ex - Ey*1i))/2;
py1 = - (Ex*M*1i)/2 + (Ey*M)/2 - (2^(1/2)*Ex*Ny)/2 - (2^(1/2)*Ey*Nx)/2 + (conj(M)*(Ey + Ex*1i))/2;
nx1 = (3*2^(1/2)*conj(Ex)*Px)/8 + (3*2^(1/2)*conj(Px)*Ex)/8 + (3*2^(1/2)*conj(Ey)*Py)/8 + (3*2^(1/2)*conj(Py)*Ey)/8;
ny1 = (2^(1/2)*conj(Ex)*Py)/8 - (2^(1/2)*conj(Ey)*Px)/8 + (2^(1/2)*conj(Px)*Ey)/8 - (2^(1/2)*conj(Py)*Ex)/8;
m1 = ((Ex - Ey*1i)*(conj(Px) - conj(Py)*1i))/8 + ((Px - Py*1i)*(conj(Ex) - conj(Ey)*1i))/8;

%     exf=fft(ex1);
%     eyf=fft(ey1);
    pxf=fft(px1);
    pyf=fft(py1);
    nxf=fft(nx1);
    nyf=fft(ny1);
    mf=fft(m1);
% ffttime = ffttime + toc;
    
%     nxf=makeReal1d(nxf);
%     nyf=makeReal1d(nyf);
    
%     epf=filter2d(epf);
%     emf=filter2d(emf);
    exf = zeros(Nr,1);
    eyf = zeros(Nr,1);
    pxf=filter1d(pxf);
    pyf=filter1d(pyf);
    nxf=filter1d(nxf);
    nyf=filter1d(nyf);
    mf=filter1d(mf);    
end