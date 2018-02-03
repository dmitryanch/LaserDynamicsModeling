clc
clear
%%
global x0 dx x
T=2e6;
tau=0.001;
Nt=ceil(T/tau);
L=2;
h=L/400;
Nr=ceil(L/h);

sigma=1.01;
gamma=2.2;
delta=0;
a=.01;
r=177.75;
% for r=[30 40]
% Npolos=10
% a=delta*(2*pi*Npolos/L)^-2
% L=2*pi*Npolos*sqrt(a/delta),    Nr=100; h=L/Nr;
l=2;
d=L;
dx=h;
% x=linspace(-d/2,d/2,Nr);
x0=(0:Nr-1)*h;
Kteor=sqrt(delta/a);Kwave=0;
t=0;
ID=round(1e5*rand*rand/rand/rand*rand/rand)/1e5;
FurT=[];
r_max=r;
sigma_min=sigma;

%%
% E0=0.003*rand(Nr,1);
% P0=0.0045*rand(Nr,1);
% D0=0.007*rand(Nr,1);

%%%%%% начальное распределение - однородный профиль
E0=sqrt(r-1-(delta/(1+sigma))^2)*ones(Nr,1)+1e-1*randn(Nr,1);
E0=sqrt(r-1-(delta/(1+sigma))^2)*ones(Nr,1);E0(1,1)=E0(1,1)+1e-2;
P0=E0*(1-1i*delta/(1+sigma));
D0=(1+(delta/(1+sigma))^2)*ones(Nr,1);

%%%%% Начальное распределение - бегущая волна
% E0=(sqrt(r-1)*exp(1i*Kteor*1.001*x0)).'+0.01*randn(Nr,1);
% P0=E0;
% D0=ones(Nr,1);

%%%% накачка - круг
% r_max=r;sigma_min=sigma;
% for m=1:Nr 
%     if (m-Nr/2-0.5)^2<=Nr^2/l^2/4
%             r(m,1)=r_max;   sigma(m,1)=sigma_min; 
%         else r(m,1)=0;sigma(m,1)=10*sigma_min;
%     end
% end
% % r=smooth(r,1);
% % sigma=smooth(sigma,1);

Pcur=P0;
Ecur=E0;
Dcur=D0;
% kx=linspace(-Nr/2,Nr/2,Nr)*2*pi/L;
kx1=(((-Nr/2):(Nr/2-1))*2*pi/L)';
kx=fftshift(kx1);
K=kx.^2;
%%
Tmax = 3e3;%+10e3;      %общее время счета
    Trecord = 1e3;      %время записи
    Tpackage = 5000;    %размер буфера для массива фазовых точек
    Out1=[];Out2=[];Out3=[];Out4=[];Out5=[];Out6=[];Out7=[];
out1 = zeros(Tpackage,1);
out2 = zeros(Tpackage,1);
out3 = zeros(Tpackage,1);
out4 = zeros(Tpackage,1);
out5 = zeros(Tpackage,1);
out6 = zeros(Tpackage,1);
out7 = zeros(Tpackage,1);
Puankare=[];
q=0;
multiplierErow=20;       %число обратное тому сколько будет браться точек из Erow
multiplierOut=10;        %число обратное тому сколько будет браться точек из Out
erec=[];
Erec=[];
outcount=1;
% OUT=zeros(Nr);
%%
tic
if t==0, PV=[];pv=[];KW=[];kw=[]; Imax=[];imax=[]; Erow=[];erow=[]; FurTT=[];FurT=[];end;
if numel(erow)>1, PV=[PV pv];pv=[];KW=[KW kw];kw=[]; Imax=[Imax imax];imax=[];Erow=[Erow erow];erow=[];FurTT=[FurTT FurT];FurT=[]; end;
    %%%% вспомогательные выражения
    exp1=exp(-gamma*tau/2);
    exp2=exp(-(1+1i*delta)*tau/2);
    exp3=exp(-(sigma+1i*a*K)*tau/2);
    exp4=sigma./(1i*a*K+sigma-1-1i*delta);
    jj=0;
    
    
while t<Tmax%*1.5
    
    oldI = abs(Ecur(1));
% tic        
    %%%% First linear part dt/2
    Dfl=exp1.*Dcur;
    
    Pfl=exp2.*Pcur;
    
    Ef=fft(Ecur);Pf=fft(Pcur).*exp4;
    Efl=ifft(exp2.*Pf+exp3.*(Ef-Pf));
    
%     %%%%%%%%%%% START Simple scheme
%     %%%% Nonliear part dt
%     Pn=Pfl+Dfl.*Efl*tau;
%     Dn=Dfl+gamma*(r-0.5*(Efl.*conj(Pfl)+Pfl.*conj(Efl)))*tau;
%     
%     %%%% Last linear part dt/2
%     Dll=exp1.*Dn;
%     Pll=exp2.*Pn;
%     Ef=fft(Efl);Pf=fft(Pn).*exp4;
%     Ell=ifft(exp2.*Pf+exp3.*(Ef-Pf));
%     
%     Ecur=Ell;   Pcur=Pll;   Dcur=Dll;
    %%%%%%%%%%%%%%%END Simple scheme

    %%%%%%%%%%%%%%%%% START Trapezional scheme
    E=Ecur;   P=Pcur;  D=Dcur;
    Ej=Ecur;
    exp5=Dcur.*Ecur;
    exp6=2*r-0.5*(conj(Ecur).*Pcur+Ecur.*conj(Pcur));
    Ef=fft(Efl);
    nmax=100;eps=1e-10; %Iteration constants
	%%% Start of Iteration Process
    for j=1:nmax
        %%%% Nonliear part dt
        Pj=Pfl+(exp5+D.*E)/2*tau;
        Dj=Dfl+gamma.*tau*[exp6-0.5*(conj(E).*P+E.*conj(P))]/2;
        Pf=fft(Pj).*exp4;
        
        %%%% Last linear part dt/2
        Dj=exp1.*Dj;
        Pj=exp2.*Pj;
        Ej=ifft(exp2.*Pf+exp3.*(Ef-Pf));

        if norm(Ej-E,inf)/norm(E,inf)<eps && norm(Pj-P,inf)/norm(P,inf)<eps && norm(Dj-D,inf)/norm(D,inf)<eps 
            E=Ej; P=Pj; D=Dj; %disp(['j=',num2str(j)]);
            break;
        else E=Ej; P=Pj; D=Dj;
        end
    end
	%%% End of Iteration Process
    if j==jj
    else jj=j;disp(['j=',num2str(j)]);
    end;
    if j==nmax
        disp(['Не сходится к ',num2str(eps),' за ',num2str(nmax),' итераций']);pause;
    end
    Ecur=E; Pcur=P; Dcur=D;
    %%%%%%%%%%%%%%%%%%%%%End Trapezional scheme
% toc 
    t=t+tau;

%       if abs(mod(round(100*t)/100,5e3))<tau/2,
%         subplot(1,2,1);plot(abs(Ecur).^2);title(['Int t=',num2str(t),'  Imax=',num2str(max(max(abs(Ecur).^2))),'  dI/Imax=',num2str((max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2)))]);xlabel('Пространственная координата x');ylabel('Уровень интенсивности I=|E|^2');
%         subplot(1,2,2);plot(real(Ecur));title([' delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);xlabel('Пространственная координата x');ylabel('Уровень реальной части поля Re(E)');
%         filename=['ID#',num2str(ID),' Int+Real 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
%         saveas(gcf,filename,'jpg');close gcf;
%     end

if(t>Tmax-Trecord)% && abs(mod(round(t/tau),10))<tau/2)
q=q+1;
out1(q) = abs(Ecur(1));
out2(q) = abs(Pcur(1));
out3(q) = Dcur(1);
out4(q) = Ecur(1)*conj(Pcur(1)) + conj(Ecur(1))*Pcur(1);
out5(q) = abs(Ecur(10));
out6(q) = (abs(Ecur(1))-oldI)/tau;
out7(q) = abs(Ecur(50));
if(q== 5000)
    Puankare =  CalcPuankare(Puankare, out1,out2,out3,out4,out5,out6,out7);
    Out1=[Out1; out1(1:multiplierOut:q)];
    Out2=[Out2; out2(1:multiplierOut:q)];
    Out3=[Out3; out3(1:multiplierOut:q)];
    Out4=[Out4; out4(1:multiplierOut:q)];
    Out5=[Out5; out5(1:multiplierOut:q)];
    Out6=[Out6; out6(1:multiplierOut:q)];
    Out7=[Out7; out7(1:multiplierOut:q)];
    q=0;
end
end;


% plot3(out1,out2,out3);drawnow;
% if(numel(out1) == 15000)
%     disp('that is all');
%     break;
% end

%     colormap(gray);
% %     set(gca,'XGrid','on');
% %     set(gca,'YGrid','on');
% %     set(gca,'ZGrid','on');
% %     %xx=linspace(0,L,Nr);
% %     %yy=linspace(0,L,Nr);
% %     %[X,Y]=meshgrid(xx,yy);
% %     %mesh(X,Y,abs(Ecur).^2);
% %     xlabel('X');
% %     ylabel('Y');
% %     zlabel('U');
% %     imagesc(atan(imag(Ecur)./real(Ecur)));
% %     imagesc(abs(fftshift(fft2(Ecur))).^2);
%     %imagesc(real(Ecur));
%     plot(1:Nr,abs(Ecur).^2);
%     title([num2str(t),'  Imax=',num2str(max(max(abs(Ecur).^2))),'  dI/Imax=',num2str((max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2)))]);
%     drawnow;
%     
%     for q=2:Nr-1
%         for w=2:Nr-1
%             if abs(Ecur(q,w))<abs(Ecur(q+1,w)) & abs(Ecur(q,w))<abs(Ecur(q,w+1)) & abs(Ecur(q,w))<abs(Ecur(q-1,w))& abs(Ecur(q,w))<abs(Ecur(q,w-1))&abs(Ecur(q,w))<0.1 %#ok<AND2>
%                Dot(q,w)=1;
%             end
%         end
%     end
% %     
%     xdz=0.05;
% %     x=x0*xdz;
%     x=xdz*linspace(-max(x0)/2,max(x0)/2,Nr);
%     Fur=dz1D(Ecur);
%     Fur = abs(fft(Ecur)).^2;
%     maxFur=max(max(Fur));for p=1:numel(Ecur), if Fur(p)==maxFur, pm=p;end,end%%%%disp(['Fur(',num2str(p),',',num2str(q),') is max']);
%     Kwave=2*pi*kx(pm);
% % % % % % % % % % %         Fur(round(numel(Ecur)/2)-10:round(numel(Ecur)/2)+10)=0;
% %     dI=(max(abs(Ecur).^2)-min(abs(Ecur).^2))/max(abs(Ecur).^2);
% %     pv=[pv dI];
%     kw=[kw Kwave];
% % % %     imax=max(abs(Ecur).^2);
%     if abs(mod(round(t/tau)*tau,5))<1e-8, FurT=[FurT Fur];  end;
    if(t>Tmax-Trecord) && abs(mod(round(t/tau),2))<tau/2
        erec=[erec Ecur(round(Nr/2),1)];
        if numel(erec)>=1e3, Erec=[Erec erec];erec=[]; 
        end
    end
    
    if abs(mod(round(t/tau),multiplierErow))<tau/2
        erow=[erow Ecur(round(Nr/2),1)];
        if numel(erow)>=1e3, Erow=[Erow erow];erow=[]; 
        end
    end
%     out(outcount)=abs( Ecur(1));
% outcount=outcount+1;
% if rem(t/tau,10) < 0.00001
% OUT=cat(2,OUT,abs(Fur));
% end
    
    
%  
% % % 
% subplot(221);
% plot(fftshift(abs(fft(Ecur)).^2));
% plot(2*pi*kx,Fur);
% % xlim([2*pi*min(x) 2*pi*max(x)]);
% % title(['t=',num2str(t),'  Imax=',num2str(max(max(abs(Ecur).^2))),'  dI/Imax=',num2str((max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))),', Kwave=',num2str(Kwave),', Kteor=',num2str(Kteor)]);
% subplot(222)
% plot(abs(Ecur).^2);
% % % subplot(223)
% % % plot(atan(imag(Ecur)./real(Ecur)));
% % title(['r=',num2str(r_max),', sigma=',num2str(sigma_min),', delta=',num2str(delta),', gamma=',num2str(gamma),', a=',num2str(a),', dx=',num2str(dx),', dt=',num2str(tau)]);
% % % % subplot(224)
% % % % plot(real(Ecur));
% drawnow;

end
toc
if numel(erow)>1, Erow=[Erow erow];erow=[];PV=[PV pv];pv=[];KW=[KW kw];kw=[]; Imax=[Imax imax];imax=[]; FurTT=[FurTT FurT];FurT=[];end;
% Ist=r-1-(delta/(1+sigma))^2;
%% Временной ряд - Локальная интенсивность
count = numel(Erow)-1;
xlimits = (numel(Erow)-count+1:numel(Erow))*tau*multiplierErow;
figure;plot(xlimits,abs(Erow(numel(Erow)-count+1:numel(Erow))).^2,'-k');
xlim([min(xlimits) max(xlimits)]);ylim([0 max(abs(Erow))^2*1.2]);
title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
xlabel('Время t');ylabel('Уровень интенсивности I=|E|^2');
filename=['ID#',num2str(ID),' Локальная интенсивность 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
% saveas(gcf,filename,'jpg');
% 
% count = 100e3;
% xlimits = (numel(Erow)-count+1:numel(Erow))*tau*multiplierErow;
% xlim([min(xlimits) max(xlimits)]);ylim([0 max(abs(Erow))^2*1.2]);
% filename=['ID#',num2str(ID),' Локальная интенсивность1 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
% saveas(gcf,filename,'jpg');
% 
% close gcf;
%%
count = 0.1*numel(Erec)-1;
xlimits = (numel(Erec)-count+1:numel(Erec))*tau*2;
figure;plot(xlimits,abs(Erec(numel(Erec)-count+1:numel(Erec))).^2,'-k');
xlim([min(xlimits) max(xlimits)]);ylim([0 max(abs(Erec))^2*1.2]);
filename=['ID#',num2str(ID),' Локальная интенсивность2 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
% saveas(gcf,filename,'jpg');
%% Профиль интенсивности
plot(abs(Ecur).^2,'-k');
ylim([0 max(abs(Ecur))^2*1.2]);
title(['r=',num2str(r_max),', sigma=',num2str(sigma_min),', delta=',num2str(delta),', gamma=',num2str(gamma),', a=',num2str(a),', dx=',num2str(dx),', dt=',num2str(tau)]);
filename=['ID#',num2str(ID),' Профиль интенсивности 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
xlabel('Координата');ylabel('Профиль интенсивности I=|E|^2');saveas(gcf,filename,'jpg');close gcf;

%% Пространственный спектр
% xdz=10;
% x0=(0:Nr-1)*h;
% x=xdz*linspace(-max(x0)/2,max(x0)/2,Nr);
% Fur=dz1D(Ecur);
% plot(2*pi*x,Fur,'-k');xlim([2*pi*min(x) 2*pi*max(x)]);
Fur=(abs(fft(Ecur)).^2);
plot(kx,Fur,'-k');
xlim([-200 200]);
filename=['ID#',num2str(ID),' Пространственные частоты 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
xlabel('Пространственная частота');ylabel('Уровень интенсивности I=|E|^2');
% saveas(gcf,filename,'jpg');close gcf;

%% Выведение спектра мощности по времени
% W1=500e3;
W1=numel(Erec);
% x0=(0:W1-1)*tau;x=x0*4;
% PTS=dz1D((Erow((numel(Erow)-W1+1):numel(Erow))).');
% figure;plot(2*pi*x,PTS,'-k');xlim([0 max(x)]);%2*pi/tau
xlimits = (1:numel(Erec((numel(Erec)-W1+1):numel(Erec))))/W1/(2*tau); %Tmax
PTS = log(abs(fft(Erec((numel(Erec)-W1+1):numel(Erec)))));
figure;plot(xlimits,PTS,'-k');
xlim([0 max(xlimits)]);
title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
xlabel('Частота');ylabel('Спектр мощности');
filename=['ID#',num2str(ID),' Спектр мощности 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
% saveas(gcf,filename,'jpg');
% xlim([0 5]);
% filename=['ID#',num2str(ID),' Спектр мощности1 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
% saveas(gcf,filename,'jpg');
% xlim([max(xlimits)-10 max(xlimits)]);
% filename=['ID#',num2str(ID),' Спектр мощности2 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
% saveas(gcf,filename,'jpg');
% close gcf;

%% Аттрактор
cortege = 1:numel(Out1);%numel(out1);
plot3(Out1(cortege),Out2(cortege),Out6(cortege),'.k');
% title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
% xlabel('abs(E)');ylabel('abs(P)');zlabel('dI');%zlabel('conj(E)*P+E*conj(P)');
filename=['ID#',num2str(ID),' Аттрактор1 delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.eps'];
% saveas(gcf,filename,'eps');close gcf;

%% Сечение Пуанкаре
section = 14;
direction = 1;
save = 1;
% Puankare=[];
% for count = 2:numel(out1);
% if count > 2 && direction * (out2(count) - section) > 0 && direction *(out2(count-1) - section)<0
%     mult = (section - out2(count-1))/(out2(count) - out2(count-1));
%     puan1 = out1(count-1) + mult * (out1(count) - out1(count-1));
%     puan3 = out3(count-1) + mult * (out3(count) - out3(count-1));
%     puan4 = out4(count-1) + mult * (out4(count) - out4(count-1));
%     puan5 = out5(count-1) + mult * (out5(count) - out5(count-1));
%     puan6 = out6(count-1) + mult * (out6(count) - out6(count-1));
%     puan7 = out7(count-1) + mult * (out7(count) - out7(count-1));
%     Puankare = [Puankare,[puan1;section;puan3;puan4;puan5;puan6;puan7]]; 
% %     plot(Puankare(1,:),Puankare(3,:),'ok'); drawnow;
% end;
% end
if(direction == 1)
    strDir = 'прямое';
else
    strDir = 'обратное';
end

figure;plot(Puankare(1,:),Puankare(3,:),'.k');
xlabel('abs(E)');ylabel('D');
title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
filename=['ID#',num2str(ID),' Сечение Пуанкаре (D) (',strDir,'-',num2str(section),') delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
if save == 1
%     saveas(gcf,filename,'jpg');close gcf;
end

figure;plot(Puankare(1,:),Puankare(4,:),'.k');
xlabel('abs(E)');ylabel('conj(E)*P+E*conj(P)');
title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
filename=['ID#',num2str(ID),' Сечение Пуанкаре (EP) (',strDir,'-',num2str(section),') delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
if save == 1
%     saveas(gcf,filename,'jpg');close gcf;
end

figure;plot(Puankare(1,:),Puankare(5,:),'.k');
xlabel('abs(E)');ylabel('abs(E1)');
title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
filename=['ID#',num2str(ID),' Сечение Пуанкаре (E1) (',strDir,'-',num2str(section),') delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
if save == 1
%     saveas(gcf,filename,'jpg');close gcf;
end
 
figure;plot(Puankare(1,:),Puankare(6,:),'.k');
xlabel('abs(E)');ylabel('dI');
title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
filename=['ID#',num2str(ID),' Сечение Пуанкаре (dI) (',strDir,'-',num2str(section),') delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
if save == 1
%     saveas(gcf,filename,'jpg');close gcf;
end
 
figure;plot(Puankare(1,:),Puankare(7,:),'.k');
xlabel('abs(E)');ylabel('abs(E2)');
title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
filename=['ID#',num2str(ID),' Сечение Пуанкаре (E2) (',strDir,'-',num2str(section),') delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
if save == 1
%     saveas(gcf,filename,'jpg');close gcf;
end
