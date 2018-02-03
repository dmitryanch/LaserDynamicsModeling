function L=matlin(E,P,D,k)
global sigma gamma delta A
L=[-sigma           A*k^2           sigma           0               0;
    -A*k^2          -sigma          0               sigma           0;
    D               0               -1              delta           real(E);
    0               D               -delta          -1              imag(E);
    -gamma*real(P)  -gamma*imag(P)  -gamma*real(E)  -gamma*imag(E)  -gamma];