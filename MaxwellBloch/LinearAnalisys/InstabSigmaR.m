clc
clear
global delta gamma a q 
tic
gamma=1;a=0.005; 
delta=-1;
for gamma=[1]
for delta=[-5]
SIGMA=[0.0001:0.0002:0.01,0.03:0.02:2,2.06:0.06:5];%;[0.00001:0.00001:0.001]
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
x5=[];
y5=[];
q=[0:0.01:0.1,0.2:0.3:30,31:1:100];
IMAG=-1*ones(length(1:1:100),length(0.05:0.05:5));q1=0;

for sigma=SIGMA %диапазон значений волнового числа
    q1=q1+1;q2=0;
    for r=1:200 %диапазон значений накачки
          q2=q2+1;qm=[];
          if r>1+delta^2/(1+sigma)^2
            Eigs=EigProgDelNegSR(sigma,r);
            i1=0;
            if all(Eigs<1e-10) %зона устойчивости
              IMAG(q2,q1)=0;
              x=[x sigma];
              y=[y r];
            elseif Eigs(1)>1e-10
                IMAG(q2,q1)=-10; %зона неустойчивости самого однородного решени€ в отсутствии модул€ции возмущений
                x3=[x3 sigma];%x2=[x2 sigma];
                y3=[y3 r];%y2=[y2 r];
                for qi=2:(numel(q)-1),
                    if Eigs(qi-1)<Eigs(qi) && Eigs(qi)>Eigs(qi+1) && Eigs(qi)>1e-10
                        i1=i1+1;qm=qi;disp(['qmax',num2str(i1),'=',num2str(q(qi)),' Increment=',num2str(Eigs(qi))]);
                    end,
                end
                if round(i1)>=1
                    IMAG(q2,q1)=-20;
                    x4=[x4 sigma];%x2=[x2 sigma];
                    y4=[y4 r];%y2=[y2 r];
                    if  all(Eigs(1:qm)>1e-10)
                        IMAG(q2,q1)=-25;
                        x5=[x5 sigma];%x2=[x2 sigma];
                        y5=[y5 r];%y2=[y2 r];
                    end
                end
            elseif all(diff(Eigs(1:2))>1e-10)
              IMAG(q2,q1)=1; %зона фазовой неустойчивости
              x1=[x1 sigma];
              y1=[y1 r];
            else %if all(Eigs(22:numel(q))<1e-10)
              IMAG(q2,q1)=2; %зона амплитудной неустойчивости
              x2=[x2 sigma];x1=[x1 sigma];
              y2=[y2 r];y1=[y1 r];
            end
          else IMAG(q2,q1)=-1;
          end
    end
end
toc
%%
figure;hold on;grid on
% border=convhull(x,y);
% patch(x(border),y(border),[0.5 0.5 0.5]);
% border2=convhull(x1,y1);
% plot(x1(border2),y1(border2),'k-.');
% x0=0:0.01:5;
% y0=1;%+(a*x2.^2-delta).^2/(1+sigma)^2;
% plot(x0,y0,'k-');

%%%%%%%% дл€ отрицательной отстройки
% X=[x1 x2 x3];Y=[y1 y2 y3];
% if numel(X)>0,border1=convhull(X,Y);
% patch(X(border1),Y(border1),[0.8 0.8 0.8]);end
% XX=[x2 x3];YY=[y2 y3];
% if numel(XX)>0, border2=convhull(XX,YY);
% patch(XX(border2),YY(border2),[0.5 0.5 0.5]);end
% if numel(x1)>0,border1=convhull(x1,y1);
% patch(x1(border1),y1(border1),[0.9 0.9 0.9]);end
if numel(x2)>0,border2=convhull(x2,y2);
patch(x2(border2),y2(border2),[0.9 0.9 0.9]);end
if numel(x3)>0,border3=convhull(x3,y3);
patch(x3(border3),y3(border3),[0.4 0.4 0.4]);end
if numel(x4)>0,border4=convhull(x4,y4);
patch(x4(border4),y4(border4),[0 0 0]);end
if numel(x5)>0,border5=convhull(x5,y5);
patch(x5(border5),y5(border5),[0.5 0 0]);end

% if numel(x1)>0,[bx11 by12]=granica(x1,y1);
% patch(bx11,by12,[0.9 0.9 0.9]);end
% if numel(x2)>0, [bx21 by22]=granica(x2,y2);
% patch(bx21,by22,[0.9 0.9 0.9]);end
% if numel(x3)>0,[bx31 by32]=granica(x3,y3);
% patch(bx31,by32,[0.4 0.4 0.4]);end
% if numel(x4)>0,[bx41 by42]=granica(x4,y4);
% patch(bx41,by42,[0 0 0]);end
% if numel(x5)>0,border5=granica(x5,y5);
% patch(x5(border5),y5(border5),[0.5 0 0]);end


hold on;plot(SIGMA,1+(delta./(1+SIGMA)).^2);
xlim([0 5]);ylim([0 200]);
title(['ќбласть нестабильности решени€ delta=',num2str(delta),' gamma=',num2str(gamma),' a=',num2str(a)]);
xlabel('Sigma');
ylabel('”ровень накачки r');
filename=['delta=',num2str(delta),' gamma=',num2str(gamma),' a=',num2str(a),' r=1..200 sigma=[',num2str(min(SIGMA)),'..',num2str(max(SIGMA)),'].jpg'];
saveas(gcf,filename,'jpg');close gcf;
figure;imagesc(IMAG);
filename=['IMAG delta=',num2str(delta),' gamma=',num2str(gamma),' a=',num2str(a),' r=1..200 sigma=[',num2str(min(SIGMA)),'..',num2str(max(SIGMA)),'].jpg'];
saveas(gcf,filename,'jpg');close gcf;
end;
end;