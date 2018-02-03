%%%% прога для учета сдвига фазы
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
for dk=[-sqrt(delta/a):0.1:-4,-3.98:0.02:4,4.1:0.1:10];% -sqrt(delta/a):0.1:10 %диапазон значений волнового числа
    om=-(2*sqrt(delta*a)*dk+a*dk.^2)/(1+sigma);
    for r=[(1-(om.^2+2*sqrt(delta*a)*om.*dk+a*om.*dk.^2)/sigma):0.1:7.9,7.91:0.01:8]   %%диапазон значений накачки
          
          Eigs=EigProgNew(dk,r);
          if all(Eigs<1e-10)
              x=[x dk];
              y=[y r];
          end
%           if all(Eigs(1:5)<1e-10) && any(Eigs>1e-10)   %%% для 1Д-неустойчивости
          if any(Eigs(1:5)>1e-10) && all(Eigs(length(q)-20:length(q))<1e-10) %%%%для 2Д-неустойчивости
              x1=[x1 dk];
              y1=[y1 r];
          end
%           q=0:0.05:0.25;
%           if  all(EigProgNew(dk,r)<-1e-10)
%               x1=[x1 dk];
%               y1=[y1 r];
%           end
    end
end
toc
figure;hold on;grid on
border2=convhull(x1,y1);
plot(x1(border2),y1(border2),'k-.');
border=convhull(x,y);
patch(x(border),y(border),[0.5 0.5 0.5]);
x2=-sqrt(delta/a):0.01:10;x2a=-(2*sqrt(delta*a)*x2+a*x2.^2)/(1+sigma);
y2=1-(x2a.^2+2*sqrt(delta*a)*x2a.*x2+a*x2a.*x2.^2)/sigma;
plot(x2,y2,'k-');

title('Область стабильности решения');
xlabel('Приращение волнового вектора k');
ylabel('Уровень накачки r');