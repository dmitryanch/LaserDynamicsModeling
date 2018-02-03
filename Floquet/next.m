function [X]=next(A,B)
global dt
A1=dt*A*B;

A2=dt*A*(B+A1/2);

A3=dt*A*(B+A2/2);

A4=dt*A*(B+A3);

X=B+(A1+2*A2+2*A3+A4)/6 ;