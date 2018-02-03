clear
clc
tic
global sigma gamma r delta A dt
sigma=5;
gamma=.1;
A=0.1;
dt=0.01;
t=0;
k=0:0.25:8;
DELTA=[0];% -0.1];%-(0:0.1:5);
R=7;
IMAG=zeros(length(R),length(DELTA));q1=0;
Ndel=20e4;Nt=5e4;
% E=sqrt(r-1.01-(delta/(1+sigma))^2);P=E*(1-1i*delta/(1+sigma));D=1+(delta/(1+sigma))^2;
a=0.1;b=0;c=0;
% E=100;P=100;D=100;
%%  решение основных уравнений
tic
for delta=DELTA %диапазон значений волнового числа
    q1=q1+1;q2=0;
    for r=R%(1+delta^2/(1+sigma)^2):1:100 %диапазон значений накачки
        q2=q2+1;
        if r>=1+delta^2/(1+sigma)^2
            t=0;
            E=zeros(1,Nt); P=zeros(1,Nt); D=zeros(1,Nt);
            E(1)=0.1;P(1)=0;D(1)=0;
            for i=1:Nt+Ndel-1
                if i<=Ndel
                    [a b c]=rk(E(1),P(1),D(1));
                    E(1)=a;P(1)=b;D(1)=c; % строка №3 для быстрого счета, без сохранения информации  
                else [a b c]=rk(E(i-Ndel),P(i-Ndel),D(i-Ndel)); 
                    E(i-Ndel+1)=a;P(i-Ndel+1)=b;D(i-Ndel+1)=c; % строка №1 для показателей ляпунова
                end
            end
%             while t<=1e3
%                 [a b c]=rk(E(numel(E)),P(numel(P)),D(numel(D)));
%                 t=t+dt;
%                 if t<=0.6e3%4e3,
%                     E=a;P=b;D=c; % строка №3 для быстрого счета, без сохранения информации  
%                 else E=[E a];P=[P b];D=[D c]; % строка №1 для показателей ляпунова
%                     if numel(E)<=1e3, else EE=[EE E];PP=[PP P];DD=[DD D];E=a;P=b;D=c;end
%                 end
%             end
%             E=EE;P=PP;D=DD;
            if max(abs(E))-min(abs(E))>0.1, IMAG(q2,q1)=2;end;
            i=1;
%             SP=zeros(5,numel(k));
            SP1=zeros(1,numel(k));
            for i=1:numel(k)
                Y=eye(5,5);Y1=[1;0;0;0;0];
                for j=1:numel(E)
                    Q=matlin(E(j),P(j),D(j),k(i)); %матрица линеаризации
                    
%                     Y=next(Q,Y); % эволюция векторов возмущения - Спектр  ляпуновских показателей
%                     [Y R]=qr(Y); %ортогонализация
%                     SP(:,i)=SP(:,i)+log(abs(diag(R)));
                    
                    Y1=next(Q,Y1); % Эволюция одного вектора возмущения - Старший ляпуновский показатель
                    Y1n=norm(Y1);
                    SP1(i)=SP1(i)+log(Y1n);
                    Y1=Y1./Y1n;
                end
            end
%             SP=SP/(numel(E)*dt);
            SP1=SP1/(numel(E)*dt); SP=SP1;
            
            %%%%% Обработка спектра старшего ляпуновского показателя
            % -1 - ниже порога
            % 0 - однородное состояние устойчиво
            % 1 - однородное состояние не устойчиво при нулевом волновом числе
            % 2 - цикл устойчив
            % 3 - цикл неустойчив при нулевом волновом числе
            % 4 - хаос при нулевых волновых числах
            % 5 - неустойчивость при ненулевых волновых числах
            % 6 - неустойчивость и при нулевых и при ненулевых волновых числах
            % 7 - несколько зон неустойчивости при ненулевых волновых числах
            k_positive=[];
            for i=1:numel(k) %%% создание списка номеров волновых чисел с положительным старшим ляпуновским показателем (СЛП)
                if SP(1,i)>0, k_positive=[k_positive i];end;
            end
            if numel(k_positive)>1 %%% подсчет количества таких зон с положительным СЛП, если они есть
                DIF=diff(k_positive);k_zone=1;k_max=1;
                for i=1:numel(k_positive)-1
                    if DIF(i)>1, k_zone=k_zone+1;end;
                end
                for i=2:numel(k_positive)-1
                    if SP(1,k_positive(i-1))<SP(1,k_positive(i)) && SP(1,k_positive(i+1))<SP(1,k_positive(i)) && round((k_positive(i+1)-k_positive(i-1)))==2
                        k_max=k_max+1;
                    end
                end
            elseif numel(k_positive)==1 && k_positive==1 && SP(1,1)<0.01, %%% если только одно волновое число с положительным СЛП
                k_zone=0;
            elseif numel(k_positive)==1 && k_positive==1 && SP(1,1)>=0.01, 
                IMAG(q2,q1)=IMAG(q2,q1)+1;
            elseif numel(k_positive)==1 && k_positive>1 && SP(1,k_positive)>=0.01
                IMAG(q2,q1)=IMAG(q2,q1)+1;
            elseif numel(k_positive)==0, 
                k_zone=0;
            end    
            if round(k_zone)==1 && all(SP(1,k_positive)<0.01) && k_max==1,
            elseif round(k_zone)==1 && k_max>1 && SP(1,1)<0.01,
                IMAG(q2,q1)=5; %%% неустойчивость при ненулевых волновых числах
            elseif round(k_zone)==1 && k_max>1 && SP(1,1)>0.01,
                IMAG(q2,q1)=6; %%% неустойчивость и при нулевых и при ненулевых волновых числах   
            elseif round(k_zone)==1 && k_positive(1)==1 && k_max==1, 
                IMAG(q2,q1)=4; %%% хаос при нулевых волновых числах;
            elseif round(k_zone)==1 && k_positive(1)>1, 
                IMAG(q2,q1)=5; %%% неустойчивость при ненулевых волновых числах
            elseif round(k_zone)>1 && k_positive(1)==1 && SP(1,1)<0.01,
                IMAG(q2,q1)=5; %%% неустойчивость при ненулевых волновых числах
            elseif round(k_zone)>1 && k_positive(1)==1 && SP(1,1)>0.01,
                IMAG(q2,q1)=6; %%% неустойчивость и при нулевых и при ненулевых волновых числах 
            elseif round(k_zone)>1,
                IMAG(q2,q1)=7; %%% несколько областей волновых чисел с положительными СЛП
            end
        else IMAG(q2,q1)=-1;
        end
    end, toc
end
toc
% figure;plot(k,SP(1,:),k,SP(2,:),k,SP(3,:),k,SP(4,:),k,SP(5,:));hold on;plot(k,zeros(numel(k)));
figure;plot(k,SP(1,:));hold on;plot(k,zeros(numel(k)));
% title(['sigma=',num2str(sigma),' gamma=',num2str(gamma),' delta=',num2str(delta),' a=',num2str(A),' r=',num2str(r)]);
% filename=['sigma=',num2str(sigma),' gamma=',num2str(gamma),' delta=',num2str(delta),' a=',num2str(A),' r=',num2str(r),' k=',num2str(min(k)),'..',num2str(max(k)),'.jpg'];
% saveas(gcf,filename,'jpg');close gcf;