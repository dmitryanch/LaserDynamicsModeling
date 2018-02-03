function [a,b,c]=rk(E,P,D)
q=0;
global dt
    a1=dt*f1(E,P) ;
    b1=dt*f2(E,P,D) ;
    c1=dt*f3(E,P,D) ;

    a2=dt*f1(E+a1/2,P+b1/2) ;
    b2=dt*f2(E+a1/2,P+b1/2,D+c1/2) ;
    c2=dt*f3(E+a1/2,P+b1/2,D+c1/2) ;

    a3=dt*f1(E+a2/2,P+b2/2) ; 
    b3=dt*f2(E+a2/2,P+b2/2,D+c2/2) ;    
    c3=dt*f3(E+a2/2,P+b2/2,D+c2/2) ;

    a4=dt*f1(E+a3,P+b3) ; 
    b4=dt*f2(E+a3,P+b3,D+c3) ;    
    c4=dt*f3(E+a3,P+b3,D+c3) ;

    
    a=E+(a1+2*a2+2*a3+a4)/6 ;
    b=P+(b1+2*b2+2*b3+b4)/6 ;
    c=D+(c1+2*c2+2*c3+c4)/6 ;
end