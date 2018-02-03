function [ef,ee]=nonlinear2dInfTetaVar(e,Ein, sigma)

    ee=ifft2(e);
    ef = filter2d(fft2(- 1i * sigma + Ein + 1i * abs(ee).^2.*ee));
    
end