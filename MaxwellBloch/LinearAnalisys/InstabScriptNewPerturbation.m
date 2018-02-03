%%%% ƒл€ возмущений нестатейного вида (см. в тетради)
clc
clear
global sigma gamma a delta q 
tic
sigma=1;
a=0.01;
delta=1;
gamma=0.01;
x=[];
y=[];
x1=[];
y1=[];
x2=[];
y2=[];
x3=[];
y3=[];
q=0:0.4:40;
K=0:0.1:20;%[0:0.1:6,6.02:0.02:14,14.1:0.1:20];
R=1:0.5:50;%,8.2:0.2:30]%[(1+(a*k^2-delta)^2/(1+sigma)^2):0.1:7.9,7.91:0.01:8] %диапазон значений накачки
IMAG=zeros(numel(R),numel(K));q1=0;
IMAG1=IMAG;

for k=K %диапазон значений волнового числа
    q1=q1+1;q2=0;
    for r=R
          q2=q2+1;
          Eigs=EigProgDelPolNew1(k,r);
          IMAG1(q2,q1)=Eigs(1);
          if all(Eigs<1e-10) %зона устойчивости
              x=[x k];
              y=[y r];              
          elseif all(Eigs(1:5)<1e-10)
              IMAG(q2,q1)=1; %зона коротковолновой неустойчивости (1D-Unstability)
              x1=[x1 k];
              y1=[y1 r];
          elseif all(Eigs(length(q)-87:length(q))<1e-10)
              IMAG(q2,q1)=2; %зона длинноволновой неустойчивости (2D-Unstability)
              x2=[x2 k];
              y2=[y2 r];
          else IMAG(q2,q1)=3;%зона длинноволновой и коротковолновой неустойчивости
              x3=[x3 k];
              y3=[y3 r];
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

title('ќбласть стабильности решени€');
xlabel('¬олновой вектор k');
ylabel('”ровень накачки r');