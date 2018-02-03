clc
clear
global sigma gamma a delta q 
tic
sigma=1;
a=0.01;
delta=1;
gamma=0.1;
x=[];
y=[];
x1=[];
y1=[];
q=0:0.1:10;

for k=[0:0.1:6,6.02:0.02:14,14.1:0.1:20] %диапазон значений волнового числа
    for r=[(1+(a*k^2-delta)^2/(1+sigma)^2):0.1:7.9,7.91:0.01:8] %диапазон значений накачки
          Eigs=EigProgExp(k,r);
          if all(Eigs<1e-10)
              x=[x k];
              y=[y r];
          end
%           if all(Eigs(1:5)<1e-10) && any(Eigs>1e-10)   %%% для 1Д-неустойчивости
          if any(Eigs(1:5)>1e-10) && all(Eigs(length(q)-20:length(q))<1e-10) %%%%для 2Д-неустойчивости
                x1=[x1 k];
                y1=[y1 r];
          end
    end
end
toc
figure;hold on;grid on
border2=convhull(x1,y1);
plot(sqrt(a)*x1(border2),y1(border2),'k-.');
border=convhull(x,y);
patch(sqrt(a)*x(border),y(border),[0.5 0.5 0.5]);
x2=0:0.01:20;
y2=1+(a*x2.^2-delta).^2/(1+sigma)^2;
plot(sqrt(a)*x2,y2,'k-');

title('Область стабильности решения');
xlabel('Волновой вектор k');
ylabel('Уровень накачки r');