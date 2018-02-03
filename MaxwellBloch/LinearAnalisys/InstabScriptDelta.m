clc
clear
global sigma gamma a q k
tic
% sigma=0.01;gamma=0.5;a=0.1;
gamma=0.001;
a=0.1;
% for gamma=[0.01 0.05 0.1 0.5 1]
for sigma=0.1%[1 0.5 0.1 0.01 0.001]%[10 5 2.5 2 1.5 1 0.5 0.1 0.01 0.001]
DELTA=-(0:0.01:1.2);
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
q=[0:0.01:0.1,0.2:0.3:30,31:100];%[0:.01:30];
% for qq=1:numel(q),qq1(qq)=-q(numel(q)+1-qq);end;q=[qq1(1:numel(q)-1) q];
IMAG=-1*ones(length(1:1:200),length(DELTA));q1=0;
%%
for delta=DELTA %диапазон значений волнового числа
    disp(['delta=',num2str(delta)]);
    q1=q1+1;q2=0;
    r=(1+delta^2/(1+sigma)^2)*1.01;
    while r<=12 %диапазон значений накачки
          q2=q2+1;qm=[];
          if r>1+delta^2/(1+sigma)^2
            Eigs=EigProgDelNeg(delta,r);
            i1=0;
            x=[x delta];
            y=[y r];
            if all(Eigs<1e-10) %зона устойчивости
              %%%
            elseif Eigs(1)>1e-10 && Eigs(1) == max(Eigs)
                x1=[x1 delta];
                y1=[y1 r];
%                 IMAG(q2,q1)=1; %Hopf instability
%                 x3=[x3 delta];x2=[x2 delta];
%                 y3=[y3 r];y2=[y2 r];
%                 for qi=2:(numel(q)-1),
%                     if Eigs(qi-1)<Eigs(qi) && Eigs(qi)>Eigs(qi+1) && Eigs(qi)>Eigs(1)%1e-10
%                         i1=i1+1;qm=qi;disp(['qmax',num2str(i1),'=',num2str(q(qi)),' Increment=',num2str(Eigs(qi))]);
%                     end,
%                 end
%                 if round(i1)>=1
% %                     IMAG(q2,q1)=-20;
%                     x4=[x4 delta];x2=[x2 delta];
%                     y4=[y4 r];y2=[y2 r];
%                     if  all(Eigs(1:qm)>1e-10)
% %                         IMAG(q2,q1)=-25;
%                         x5=[x5 delta];x2=[x2 delta];
%                         y5=[y5 r];y2=[y2 r];
%                     end
%                 end
%             elseif all(diff(Eigs(1:2))>1e-10)
%               IMAG(q2,q1)=1; %зона фазовой неустойчивости
%               x1=[x1 delta];
%               y1=[y1 r];
            elseif Eigs(1)<1e-10
%               IMAG(q2,q1)=2; % wave instability
              x2=[x2 delta];%x1=[x1 delta];
              y2=[y2 r];%y1=[y1 r];
            else
                x3=[x3 delta];x2=[x2 delta];x1=[x1 delta];
                y3=[y3 r];y2=[y2 r];y1=[y1 r];
            end
%           else IMAG(q2,q1)=-1;
          end
          r=r+0.0125;
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

%%%%%%%% для отрицательной отстройки
% X=[x1 x2 x3];Y=[y1 y2 y3];
% if numel(X)>0,border1=convhull(X,Y);
% patch(X(border1),Y(border1),[0.8 0.8 0.8]);end
% XX=[x2 x3];YY=[y2 y3];
% if numel(XX)>0, border2=convhull(XX,YY);
% patch(XX(border2),YY(border2),[0.9 0.9 0.9]);end
% if numel(x3)>0,border3=convhull(x3,y3);
% patch(x3(border3),y3(border3),[0.4 0.4 0.4]);end
% if numel(x4)>0,border4=convhull(x4,y4);
% patch(x4(border4),y4(border4),[0 0 0]);end

if numel(x)>0,[bx1 by2]=granica(x,y);
patch(bx1,by2,[0.9 0.9 0.9],'LineWidth',1);end
if numel(x1)>0,[bx11 by12]=granica(x1,y1);
patch(bx11,by12,[0.7 0.7 0.7],'LineWidth',1);end
if numel(x2)>0, [bx21 by22]=granica(x2,y2);
patch(bx21,by22,[0.6 0.6 0.6],'LineWidth',1);end
if numel(x3)>0,[bx31 by32]=granica(x3,y3);
patch(bx31,by32,[0.4 0.4 0.4],'LineWidth',1);end
if numel(x4)>0,[bx41 by42]=granica(x4,y4);
patch(bx41,by42,[0 0 0],'LineWidth',1);end
if numel(x5)>0,[bx51 by52]=granica(x5,y5);
patch(bx51,by52,[0.5 0 0],'LineWidth',1);end

hold on;plot(-(0.0:0.1:10),1+(-(0.05:0.1:10)./(1+sigma)).^2,'black','LineWidth',2);
xlim([-5 max(DELTA)]);ylim([0 100]);
% title(['Область нестабильности решения sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a)]);
set(gca,'FontSize',16);
xlabel('Отстройка частоты \delta','FontSize',18,'FontWeight','bold');
ylabel('Уровень накачки r','FontSize',18,'FontWeight','bold');
filename=['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=0..200 delta=',num2str(min(DELTA)),'..',num2str(max(DELTA))];
% saveas(gcf,[filename,'.jpg'],'jpg');%close gcf;
% save([filename,'.mat'])
% figure;imagesc(IMAG);colorbar('vert');
% filename=['IMAG sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=0..200 delta=',num2str(min(DELTA)),'..',num2str(max(DELTA)),'.jpg'];
% saveas(gcf,filename,'jpg');close gcf;

end;
% end;

%% Диаграмма всех собственных чисел для всех q
% delta=-1;r=30;figure;[X Y]=EigProgDelNeg(delta,r);plot(q,X);title(['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)]);Y(:,(numel(q)-1)/2+1)
% % k=11.3;figure;[X Y]=EigProgDelNegRandK(delta,r);plot(q,X);title(['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)]);Y(:,(numel(q)-1)/2+1)
% %%% Нахождение всех максимумов наибольшего собственного числа
% % % i1=0;for i=2:(numel(q)-1),if X(5,i-1)<X(5,i) && X(5,i)>X(5,i+1),i1=i1+1;disp(['qmax',num2str(i1),'=',num2str(q(i)),' i=',num2str(i),' Increment=',num2str(X(5,i))]);end,end
% disp(['k=',num2str(k)]);
% i1=0;for i=2:(numel(q)-1),if X(i-1)<X(i) && X(i)>X(i+1)&& X(i)>1e-10,i1=i1+1;disp([' qmax',num2str(i1),'=',num2str(q(i)),' ===> Knew=',num2str(k+q(i)),' i=',num2str(i),' Increment=',num2str(Y(i))]);end,end
