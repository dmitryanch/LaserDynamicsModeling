function [ef,ee]=nonlinear2dInf(e,Ein)

    ee=ifft2(e);
    ef = filter2d(fft2(Ein + 1i * abs(ee).^2.*ee));
    
end