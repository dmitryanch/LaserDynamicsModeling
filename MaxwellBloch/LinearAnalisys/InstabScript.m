clc
clear
global sigma gamma a delta q 
tic
sigma=2.2;
a=0.01;
delta=1;
gamma=1;
x=[];
y=[];
x1=[];
y1=[];
x2=[];
y2=[];
x3=[];
y3=[];
x4=[];
y4=[];
q=0:5:500;
K=0:0.1:20;%[0:0.1:6,6.02:0.02:14,14.1:0.1:20];
R=1:0.05:10;%,8.2:0.2:30]%[(1+(a*k^2-delta)^2/(1+sigma)^2):0.1:7.9,7.91:0.01:8] %диапазон значений накачки
IMAG=zeros(numel(R),numel(K));q1=0;
IMAG1=IMAG;

for k=K %диапазон значений волнового числа
    q1=q1+1;q2=0;
    for r=(1+(a*k^2-delta).^2/(1+sigma)^2):0.05:5
          q2=q2+1;
          Eigs=EigProg(k,r);
          IMAG1(q2,q1)=Eigs(1);
          if all(Eigs<1e-10) %зона устойчивости
              x=[x k];
              y=[y r];              
          elseif Eigs(1)>1e-10 % временная неустойчивость решения
              IMAG(q2,q1)=-2; %зона неустойчивости самого однородного решения в отсутствии модуляции возмущений
                x4=[x4 k];
                y4=[y4 r];
          elseif all(diff(Eigs(1:2))<1e-10)    % 1D- inStab - /амплитудная неустойчивость/, 2D-inStab - /фазовая неустойчивость/
              IMAG(q2,q1)=1; 
              x1=[x1 k];
              y1=[y1 r];
          else                                 % 1D- inStab - /фазовая неустойчивость/, 2D-inStab - /амплитудная неустойчивость/
              IMAG(q2,q1)=2;
              x2=[x2 k];
              y2=[y2 r];
          end
    end
end
toc
%%
figure;hold on;grid on
border2=convhull(x1,y1);
plot(sqrt(a)*x1(border2),y1(border2),'k-.');
border=convhull(x,y);
patch(sqrt(a)*x(border),y(border),[0.5 0.5 0.5]);
x0=0:0.01:20;
y0=1+(a*x0.^2-delta).^2/(1+sigma)^2;
plot(sqrt(a)*x0,y0,'k-');

title('Область стабильности решения');
xlabel('Волновой вектор k');
ylabel('Уровень накачки r');