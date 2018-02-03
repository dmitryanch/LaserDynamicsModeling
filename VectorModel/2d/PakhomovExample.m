clc 
clear

D=1;
a=0.01;
gamma=1;
sigma=0.1;
delta=1;
omega=delta;
k=sqrt(delta/a);

r=1:0.5:25;
q=0:0.5:25;

x=[];
y=[];
Res=-1*ones(numel(r),numel(q));
for i=1:numel(q)
    for j=1:numel(r)
        E=sqrt(r(j)-1);
        P=E;
        J=[+1i*omega+sigma+1i*a*(k+q(i))^2, 0, -sigma, 0, 0;
            0, -1i*omega+sigma-1i*a*(k-q(i))^2, 0, -sigma, 0;
            -D, 0, 1i*omega+(1+1i*delta),0,-E;
            0, -D, 0, -1i*omega+(1-1i*delta), -E;
            1/2*conj(P), 1/2*P, 1/2*conj(E), 1/2*E, gamma];
%         e=eig(J);
        lyam=sym('lyam');
        J=J+lyam*eye(5);
        e=double(solve(det(J),lyam));
        if max(real(e))>1e-8
            Res(j,i)=1;
            x=[x q(i)];
            y=[y r(j)];
            plot(x,y);
            drawnow;
            find=1;
%             break;
        else
            Res(j,i)=0;
        end
    end
end
