function dispMode(omega,a,sigma,b,c,gA,gP,r)

kx0=real((omega - gP)^(1/2))/a^(1/2);
ky0=real((gP + omega)^(1/2))/a^(1/2);
disp(['Kx=',num2str(kx0),' Ky=',num2str(ky0)])
rthx=(gA + sigma)/sigma + (a*kx0^2 + gP - omega)^2/(sigma*(gA + sigma + 1)) - (a*kx0^2 + gP - omega)^2/(sigma*(gA + sigma + 1)^2);
rthy=(- a*ky0^2 + gP + omega)^2/(sigma*(sigma - gA + 1)) - (gA - sigma)/sigma - (- a*ky0^2 + gP + omega)^2/(sigma*(sigma - gA + 1)^2);
disp(['Rthx=',num2str(rthx),' Rthy=',num2str(rthy)])
wx0 = omega + (a*kx0^2 + gP - omega)/(gA + sigma + 1);
wy0 = omega - (- a*ky0^2 + gP + omega)/(sigma - gA + 1);
disp(['Wx=',num2str(wx0),' Wy=',num2str(wy0)])
Ix = (4*b*c*r*sigma)/((gA + sigma)*(b + 3*c)) - (4*b*c*(a*kx0^2 + gP - omega)^2)/((b + 3*c)*(gA + sigma + 1)^2) - (4*b*c)/(b + 3*c);
Iy = - (4*b*c)/(b + 3*c) - (4*b*c*(- a*ky0^2 + gP + omega)^2)/((b + 3*c)*(sigma - gA + 1)^2) - (4*b*c*r*sigma)/((b + 3*c)*(gA - sigma));
disp(['Ix=',num2str(Ix),' Iy=',num2str(Iy)])