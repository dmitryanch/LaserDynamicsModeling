clc
clear
%%
global x0 dx x
T=8e3;
t=0;
tau=0.02;
Nt=ceil(T/tau);
sigma=0.1;
gamma=0.005;
delta=5;
a=.0001;
r=4;
diag=1;
Npolos=0;
L=1;
l=0.75;
NrecSize = 0;   %% parameter for cutting informational part of profile
if delta > 0 && Npolos>0
    if diag == 1
        L=2*pi/sqrt(delta/a/2)*Npolos/l;
    else 
        L=2*pi/sqrt(delta/a)*Npolos/l;
    end
else
    L=L/l;
end
h=L/256;
Nr=ceil(L/h);
stopfilename='stop';
kfil=int32((Nr-1)/3)+1;
% int Nr = 256, Nt=1;
%     double t = 0.0, tau = 0.01, L = 5.0, h = L/Nr;		
%     double tau_2 = tau / 2;
%     double gamma =0.1;
%     double sigma =10;
%     double delta = -4.5;
%     double a =0.01;
%     double r = 100.0;

r_max=r;sigma_min=sigma;

% Npolos=8;disp(['количество полос на апертуре: ',num2str(Npolos)]);
% a=delta*(2*pi*Npolos/L)^-2;

d=L;
dx=h;
% x0=linspace(-d/2,d/2,Nr);
x0=(0:Nr-1)*h;
[X Y]=meshgrid(x0);
Kt=sqrt(delta/a);Kwave=0;

%%% накачка - квадрат
r=-1*ones(Nr,Nr);sigma=sigma_min*99*ones(Nr,Nr);
sigma(round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)),round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)))=sigma_min;
r(round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)),round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)))=r_max;
r=smooth(r,1);
sigma=smooth(sigma,1);

%%%% накачка - круг
% r_max=r;sigma_min=sigma;
% for m=1:Nr 
%     for n=1:Nr
%         if (m-Nr/2-0.5)^2+(n-Nr/2-0.5)^2<=Nr^2/l^2/4
%             r(m,n)=r_max;   sigma(m,n)=sigma_min; 
%         else r(m,n)=0;sigma(m,n)=5*sigma_min;
%         end
%     end
% end
% r=smooth(r,1);
% sigma=smooth(sigma,1);
%%
% E0=0.3*rand(Nr,Nr);
% P0=0.4*rand(Nr,Nr);
% D0=0.5*rand(Nr,Nr);
% E0=smooth(0.001*randn(Nr,Nr),1);
% P0=smooth(0.001*randn(Nr,Nr),1);
% D0=smooth(0.001*randn(Nr,Nr),1);

%%%%%% начальное распределение - однородный профиль
E0=(sqrt(r_max-1-(delta./(1+sigma_min)).^2)/sqrt(2)+1i*sqrt(r_max-1-(delta./(1+sigma_min)).^2)/sqrt(2)).*r/r_max;%+0.001*randn(Nr,Nr);
P0=E0*(1-1i*delta./(1+sigma_min));
D0=(1+(delta./(1+sigma_min)).^2).*r/r_max;
% E0(1,1)=E0(1,1)+1e-10;

%%%%%% начальное распределение - квадратна€ решетка вихрей
% A1=sqrt((r_max-1)/5);A2=A1;A3=A1;A4=A1;
% f1=pi;f2=pi;f3=pi;f4=0;
% if diag == 1
%     kx=Kt/sqrt(2);ky=Kt/sqrt(2);E0=A1*exp(1i*(-kx*X-ky*Y+f1)) + A2*exp(1i*(kx*X+ky*Y+f2)) + A3*exp(1i*(kx*X-ky*Y+f3)) + A4*exp(1i*(-kx*X+ky*Y+f4));
% else
%     kx=Kt;ky=Kt;E0=A1*exp(1i*(-kx*X+f1)) + A2*exp(1i*(kx*X+f2)) + A3*exp(1i*(-ky*Y+f3)) + A4*exp(1i*(ky*Y+f4));
% end
% P0=E0;
% D0=ones(Nr,Nr).*(r/5+4/5)-A1^2*exp(1i*(-2*kx*X)) - A2^2*exp(1i*(2*kx*X)) + A3^2*exp(1i*(-2*ky*Y)) + A4^2*exp(1i*(2*ky*Y));

Pcur=P0;
Ecur=E0;
Dcur=D0;
% kx=linspace(-Nr/2,Nr/2,Nr)*2*pi/L;
kx=((-Nr/2):(Nr/2-1))*2*pi/L;
% [X Y]=meshgrid(kx);
% K=fftshift(X.^2+Y.^2);
kx=fftshift(kx);
ky=kx;
% ky=fftshift(ky);
K=zeros(Nr,Nr);
for i=1:Nr
    for j=1:Nr
        K(i,j)=kx(i)^2+ky(j)^2;
    end
end
nstep=int32(0);
nrec=int32(0);
erow_incr=10;
ID=int16(0);
datafilename=['INT delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.dat'];
while exist(datafilename, 'file')==2,
    ID=ID+1;
    datafilename=['INT delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.dat'];
end
nl=1000;
%%
tic
% Dot=zeros(Nr,Nr);
if t==0, 
    erow=zeros(nl,1); %Erow=[];DI=[];dI=[]; Imax=[];imax=[]; 
end;
% if numel(erow)>1, Erow=[Erow erow];erow=[]; DI=[DI dI];dI=[]; Imax=[Imax imax];imax=[]; end;
    %%%% вспомогательные выражени€
    exp1=exp(-gamma*tau/2);
    exp2=exp(-(1+1i*delta)*tau/2);
    exp3=exp(-1i*a*K*tau/2);
    
    jj=0;
    F=[];
    

while round(t/tau)*tau<20e4*tau%6e4
    
    nstep=nstep+int32(1);
% tic        
    %%%% First linear part dt/2
    Dfl=exp1.*Dcur;
    Pfl=exp2.*Pcur;
    Efl=ifft2(exp3.*fft2(Ecur));
    
    %%%%%%%%%%%%%%%%% START Trapezional scheme
    E=Ecur;   P=Pcur;  D=Dcur;
    Ej=Ecur;
    exp4=sigma.*(Pcur-Ecur);
    exp5=Dcur.*Ecur;
    exp6=2*r-0.5*(conj(Ecur).*Pcur+Ecur.*conj(Pcur));
%     Ef=fft2(Efl);
    nmax=100;eps=1e-10; % Iteration Constants
	%%%% Start of Iteration Process
    for j=1:nmax
        %%%% Nonlinear part dt
        Ej=Efl+(exp4+sigma.*(P-E))/2*tau;
        Pj=Pfl+(exp5+D.*E)/2*tau;
        Dj=Dfl+gamma.*[exp6-0.5*(conj(E).*P+E.*conj(P))]/2*tau;
        
        %%%% Last linear part dt/2
        Dj=exp1.*Dj;
        Pj=exp2.*Pj;
        Ej=ifft2(exp3.*fft2(Ej));
        
        
        if norm(Ej-E,inf)/norm(E,inf)<eps && norm(Pj-P,inf)/norm(P,inf)<eps && norm(Dj-D,inf)/norm(D,inf)<eps 
            E=Ej; P=Pj; D=Dj;
            %disp(['j=',num2str(j)]);
            break;
        else E=Ej; P=Pj; D=Dj;
        end
    end
	%%%% End of Iteration Process
    if j==jj
    else jj=j;disp(['j=',num2str(j)]);
    end;
    if j==nmax
        disp(['Ќе сходитс€ к ',num2str(eps),' за ',num2str(nmax),' итераций']);%pause;
    end
    Ecur=E; Pcur=P; Dcur=D;
    %%%%%%%%%%%%%%%%%%%%%End Trapezional scheme
% toc 
    t=t+tau;
    if abs(mod(nstep,int32(erow_incr)))<tau/2,
%         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
%         disp(t);
        nrec=nrec+1;
        erow(nrec)=Ecur(Nr/2,Nr/2);%Ecur(round(Nr/2),round(Nr/2))];
        
         if nrec==nl, 
    %         Erow=[Erow erow];erow=[];%DI=[DI dI];dI=[]; % Imax=[Imax imax];imax=[];
    %         disp('write center');
            disp(['t= ',num2str(t),' at ',num2str(toc),' s']);
            fid = fopen(datafilename, 'a');
            fwrite(fid, abs(erow).^2, 'double');
            fclose(fid); 
            nrec=int32(0);   % запись матрицы в файл (40 байт)
            if exist(stopfilename,'file')==2,
                break;
            end
        end
    end
   
%     if abs(mod(round(100*t)/100,3e2))<tau/2,
%         imagesc(abs(Ecur).^2);title(['Int t=',num2str(t),'  Imax=',num2str(max(max(abs(Ecur).^2))),'  dI/Imax=',num2str((max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))),' delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a)]);
%         filename=['Int2D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' ID#',num2str(ID),'.jpg'];
%         saveas(gcf,filename,'jpg');close gcf;
%     end
    
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
%     imagesc(atan(imag(Ecur)./real(Ecur)));
% %     imagesc(abs(fftshift(fft2(Ecur))).^2);
%     %imagesc(real(Ecur));
%     imagesc(abs(Ecur).^2);colorbar('vert');
% % % contour(real(Ecur),1);hold on;contour(imag(Ecur),1);hold off;
%     title([num2str(t),'  Imax=',num2str(max(max(abs(Ecur).^2))),'  dI/Imax=',num2str(dI(numel(dI)))]);
%     colormap(gray);
%     drawnow;
%     
%     for q=2:Nr-1
%         for w=2:Nr-1
% %             if abs(Ecur(q,w))< & abs(Ecur(q,w))<abs(Ecur(q,w+1)) & abs(Ecur(q,w))<abs(Ecur(q-1,w))& abs(Ecur(q,w))<abs(Ecur(q,w-1))&abs(Ecur(q,w))<0.1 %#ok<AND2>
%             if abs(Ecur(q,w))==min(min(abs(Ecur(q-1:q+1,w-1:w+1)))) && abs(Ecur(q,w))<0.2
%                 Dot(q,w)=1;
%             end
%         end
%     end
%     
%     xdz=0.8;
% % %     x=x0*xdz;
%     x=xdz*linspace(-max(x0)/2,max(x0)/2,Nr);
%     Fur=dz(Ecur);
%     maxFur=max(max(Fur));for p=1:size(Ecur,1), for q=1:size(Ecur,2), if Fur(p,q)==maxFur, pm=p;qm=q;break;end,end,end%%%%disp(['Fur(',num2str(p),',',num2str(q),') is max']);
%     Kwave=2*pi*sqrt((x(pm))^2+(x(qm))^2);
% % % % %         Fur(round(size(Ecur,1)/2),round(size(Ecur,1)/2))=0;
% % 
% subplot(221);
% % imagesc(fftshift(abs(fft2(Ecur))^2));
% imagesc(Fur);
% title(['t=',num2str(t),'  Imax=',num2str(max(max(abs(Ecur).^2))),'  dI/Imax=',num2str((max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))),', Kwave=',num2str(Kwave),', Kteor=',num2str(Kt)]);
% subplot(222)
% imagesc(abs(Ecur).^2);colorbar('vert');
% % mesh(abs(Ecur).^2);
% subplot(223)
% imagesc(atan(imag(Ecur)./real(Ecur)));
% title(['r=',num2str(r),', sigma=',num2str(sigma),', delta=',num2str(delta),', gamma=',num2str(gamma),', a=',num2str(a),', dx=',num2str(dx),', dt=',num2str(tau)]);
% subplot(224)
% % imagesc(real(Ecur));
% % imagesc(-Dot);Dot=Dot*0.99;%title(['t=',num2str(t),'  Imax=',num2str(max(max(abs(Ecur).^2))),'  dI/Imax=',num2str((max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))),', Kwave=',num2str(Kwave),', Kteor=',num2str(Kt)]);
% colormap(gray);
% drawnow;

%%%%% Video part
% if abs(mod(round(t/tau)*tau,30))<1e-8
%     imagesc(atan(imag(Ecur)./real(Ecur)));
%     colormap(gray);
%     filename=['faz',num2str(round(t/tau)),'.jpg'];
%     saveas(gcf,filename,'jpg');
%     
%     imagesc(abs(Ecur).^2);
%     colormap(gray);
%     filename=['bz',num2str(round(t/tau)),'.jpg'];
%         saveas(gcf,filename,'jpg');
% 
%     contour(real(Ecur),1e-15);hold on;contour(imag(Ecur),-1e-15);hold off;
%     colormap(winter);
%     filename=['con',num2str(round(t/tau)),'.jpg'];
%     saveas(gcf,filename,'jpg');
%     
%     xdz=40;
% % %     x=x0*xdz;
%     x=xdz*linspace(-max(x0)/2,max(x0)/2,Nr);
%     Fur=dz(Ecur);
%     imagesc(Fur);
%     colormap(gray);
%     filename=['dz',num2str(round(t/tau)),'.jpg'];
%     saveas(gcf,filename,'jpg');
%     for p=1:size(Ecur,1), for q=1:size(Ecur,2), if 2*pi*sqrt(x(p)^2+x(q)^2)<5, Fur(p,q)=0;end,end,end
%      imagesc(Fur);
%     colormap(gray);
%     filename=['1dz',num2str(round(t/tau)),'.jpg'];
%     saveas(gcf,filename,'jpg');
    
%     imagesc(-Dot);Dot=Dot*0.99;
%     colormap(gray);
%     filename=['dots',num2str(round(t/tau)),'.jpg'];
%     saveas(gcf,filename,'jpg');
% end
end
toc;
%%
if nrec>0, 
%     Erow=[Erow erow];erow=[];  DI=[DI dI];dI=[]; Imax=[Imax imax];imax=[]; 
%     disp('write end');
    fid = fopen(datafilename, 'a');    
    fwrite(fid, abs(erow(1:nrec)).^2, 'double');
    fclose(fid); 
%     disp(nrec);
    nrec=int32(0);
    erow=zeros(nl,1);
end;

%%
INT=[];
fid = fopen(datafilename, 'rb');
INT = fread(fid, 'double'); 
fclose(fid); 
% Ist=r_max-1-(delta/(1+sigma_min))^2;
figure;plot((1:numel(INT))*erow_incr*tau,INT);
title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
filename=['Int(t) 2D delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.jpg'];
% xlabel('¬рем€ t');ylabel('”ровень интенсивности I=|E|^2');saveas(gcf,filename,'jpg');close gcf;

%% ¬ыведение спектра мощности по времени
% W1=2000;
% x0=(0:W1-1)*tau;x=x0*0.1;
% PTS=dz1D((Erow((numel(Erow)-W1+1):numel(Erow))).');
% PTS=dz1D((Erow(6.45e5:6.47e5-1)).');
% figure;plot(2*pi*x,PTS);xlim([0 50]);%2*pi/tau
% title(['delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L)]);
% filename=['ID#',num2str(ID),' Power time-Spectra 1D delta=',num2str(delta),' r=',num2str(r),' sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' t=',num2str(t),'.jpg'];
% xlabel('¬рем€ t');ylabel('”ровень интенсивности I=|E|^2');saveas(gcf,filename,'jpg');close gcf;
%%
x0=(0:Nr-1)*h;
gray=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';
figure;imagesc(x0(Nr*NrecSize+1:Nr*(1-NrecSize)),x0(Nr*NrecSize+1:Nr*(1-NrecSize)),abs(Ecur(Nr*NrecSize+1:Nr*(1-NrecSize),Nr*NrecSize+1:Nr*(1-NrecSize))).^2);
colorbar('vert');     
colormap(gray);
% filename=['bz I ',num2str(round(t/tau)),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% saveas(gcf,filename,'png');
% filename=['bz I ',num2str(round(t/tau)),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.eps'];
% saveas(gcf,filename,'eps');
% close gcf;
% xdz=100;    x=xdz*linspace(-max(x0)/2,max(x0)/2,Nr);    Fur=dz(Ecur); figure;   imagesc(2*pi*x(Nr*3/16:Nr*13/16),2*pi*x(Nr*3/16:Nr*13/16),Fur(Nr*3/16:Nr*13/16,Nr*3/16:Nr*13/16)); colorbar('vert'); colormap(gray);
% filename=['dz I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% saveas(gcf,filename,'png');
% filename=['dz I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.eps'];
% saveas(gcf,filename,'eps');
% maxFur=max(max(Fur));for p=1:size(Ecur,1), for q=1:size(Ecur,2), if Fur(p,q)==maxFur, pm=p;qm=q;break;end,end,end%%%%disp(['Fur(',num2str(p),',',num2str(q),') is max']);
% Kwave=2*pi*sqrt((x(pm))^2+(x(qm))^2)
% close gcf;