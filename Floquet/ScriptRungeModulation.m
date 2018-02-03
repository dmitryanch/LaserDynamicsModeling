clear
clc
tic
global sigma gamma r delta A dt
sigma=0.1;
gamma=.001;
delta=-0.25;
r0=5;
A=0.0001;
dt=0.5;
t=0;
m=0.25;
wRel=sqrt(2*gamma*sigma*(r0-1-delta^2)/(1+delta^2));
w=1*wRel;

% E=0.1;P=0;D=0;
E=sqrt(r0-1.01-(delta/(1+sigma))^2);P=E*(1-1i*delta/(1+sigma));D=1+(delta/(1+sigma))^2;
% E=sin(2)+1i*cos(2);P=sin(4)+1i*cos(4);D=tan(5);
%%  решение основных уравнений
q=1;EE=[];PP=[];DD=[];
tic
while t<=1.5e6*dt
    r=r0*(1+m*sin(w*t));
    [a b c]=rk(E(numel(E)),P(numel(P)),D(numel(D)));
    t=t+dt;q=q+1;
    if t<=1e6*dt%1e3%4e3,
        E=a;P=b;D=c; % строка №3 для быстрого счета, без сохранения информации  
%         E=[E a];P=b;D=c; % строка №2 для временного ряда интенсивности
    else E=[E a];P=[P b];D=[D c]; % строка №1 для показателей ляпунова
        if numel(E)<=1e3, else EE=[EE E];PP=[PP P];DD=[DD D];E=a;P=b;D=c;end
    end
%     plot((0:numel(E)-1)*dt,abs(E).^2);drawnow;
end
if numel(E)>1, EE=[EE E];PP=[PP P];DD=[DD D]; E=a;P=b;D=c;end;
figure;plot((0:numel(EE)-1)*dt,abs(EE).^2);drawnow;
toc
%% ячейка с показателями ляпунова (ветора возмущения)
E=EE;P=PP;D=DD;
% E=EE(numel(EE)-8e5:numel(EE));P=PP(numel(PP)-8e5:numel(PP));D=DD(numel(DD)-8e5:numel(DD));
% E=E(round(numel(E)/2):numel(E));P=P(round(numel(P)/2):numel(P));D=D(round(numel(D)/2):numel(D));
k=0:0.1:20;
i=1;SP=zeros(5,numel(k));
tic
    for i=1:numel(k)
        Y=eye(5,5);
        for j=1:numel(E)
            Q=matlin(E(j),P(j),D(j),k(i)); %матрица линеаризации
            Y=next(Q,Y); % эволюция векторов возмущения
            [Y R]=qr(Y); %ортогонализация
            SP(:,i)=SP(:,i)+log(abs(diag(R)));
        end
        disp([num2str(i/numel(k)*100),'%']);        
    end
toc
SP1=SP/(numel(E)*dt);

%%
figure;plot(k,SP1(1,:),'black','LineWidth',1);hold on;grid on;
plot(k,SP1(2,:),'black','LineWidth',1)
plot(k,SP1(3,:),'black','LineWidth',1)
plot(k,SP1(4,:),'black','LineWidth',1)
plot(k,SP1(5,:),'black','LineWidth',1);
hold on;plot(k,zeros(numel(k)),'black');
% ylim([-0.5 0.1])
% title(['sigma=',num2str(sigma),' gamma=',num2str(gamma),' delta=',num2str(delta),' a=',num2str(A),' r=',num2str(r)]);
filename=['sigma=',num2str(sigma),' gamma=',num2str(gamma),' delta=',num2str(delta),' a=',num2str(A),' r=',num2str(r),' k=',num2str(min(k)),'..',num2str(max(k)),'.jpg'];
set(gca,'FontSize',16);
xlabel('q','FontSize',18,'FontWeight','bold');
ylabel('\lambda','FontSize',18,'FontWeight','bold');
% saveas(gcf,filename,'jpg');close gcf;