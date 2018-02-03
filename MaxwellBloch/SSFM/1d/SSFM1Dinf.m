%%%% Simple script for MB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are can be set as on infinite or on finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON delta_aperture_ON sigma_aperture_ON
T=4e5;
t=0;
tau=0.1;
Nt=ceil(T/tau);
sigma=0.05;
gamma=0.0005;
delta=-1;
a=.01;
r=10;

%%%% indicates parameter profiling
delta_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;

%%%%% Transverse shape of active medium
apertureType = 'square';
% apertureType = 'circle';

diag=1;    %% use diagonal simmetry of square vortex latticces or 
Npolos=0;
L=50;
l=0.75;         %% active part of grid
NrecSize = 0;   %% parameter for cutting informational part of profile
if delta > 0 && Npolos > 0
    if diag == 1
        L=2*pi/sqrt(delta/a/2)*Npolos/l;
    else 
        L=2*pi/sqrt(delta/a)*Npolos/l;
    end
else
    if delta_aperture_ON==1 || sigma_aperture_ON==1 || r_aperture_ON==1
        L=L/l;
    end
end
h=L/1024;
Nr=ceil(L/h);
stopfilename='stop0';        %%% file name that will stop calculating
kfil=int32((Nr-1)/7*3)+1;     %%% parameter for dealising

kx=((-Nr/2):(Nr/2-1))*2*pi/L;
kx=fftshift(kx);
K=kx.^2;

nstep=int32(0);
nrec=int32(0);
erow_incr=10;
nl=1000;

nw=8192;w=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);

disp(['Kx_max=',num2str(max(kx))]);
disp(['dk=',num2str(kx(2))]);
disp(['W_max=',num2str(max(w))]);
disp(['dW=',num2str(2*pi/(tau*nw*erow_incr))]);
disp(['OMEGA_0_teor=',num2str(sigma*delta/(1+sigma))]);

r_max=r;sigma_min=sigma;

% Npolos=8;disp(['количество полос на апертуре: ',num2str(Npolos)]);
% a=delta*(2*pi*Npolos/L)^-2;

d=L;
dx=h;
% x0=linspace(-d/2,d/2,Nr);
x0=(0:Nr-1)*h;
[X Y]=meshgrid(x0);
% Kt=sqrt(delta/a);Kwave=0;
Nt=100e6;
Kt=Nt/nl/erow_incr;
KK=zeros(numel(kx),Kt);
nk=1;  
%%
% E0=0.3*rand(Nr,Nr);
% P0=0.4*rand(Nr,Nr);
% D0=0.5*rand(Nr,Nr);
% E0=smooth(0.001*randn(Nr,Nr),1);
% P0=smooth(0.001*randn(Nr,Nr),1);
% D0=smooth(0.001*randn(Nr,Nr),1);

% if delta >= 0
%     %%%%%% начальное распределение - квадратная решетка вихрей
%     A1=sqrt((r_max-1)/5);A2=A1;A3=A1;A4=A1;
%     f1=pi;f2=pi;f3=pi;f4=0;
%     if diag == 1
%         kx=Kt/sqrt(2);ky=Kt/sqrt(2);E0=A1*exp(1i*(-kx*X-ky*Y+f1)) + A2*exp(1i*(kx*X+ky*Y+f2)) + A3*exp(1i*(kx*X-ky*Y+f3)) + A4*exp(1i*(-kx*X+ky*Y+f4));
%         D0=ones(Nr,Nr).*(r/5+4/5)-A1^2*exp(1i*(-2*kx*X-2*ky*Y)) - A2^2*exp(1i*(2*kx*X+2*ky*Y)) + A3^2*exp(1i*(2*kx*X-2*ky*Y)) + A4^2*exp(1i*(-2*kx*X+2*ky*Y));
%     else
%         kx=Kt;ky=Kt;E0=A1*exp(1i*(-kx*X+f1)) + A2*exp(1i*(kx*X+f2)) + A3*exp(1i*(-ky*Y+f3)) + A4*exp(1i*(ky*Y+f4));
%         D0=ones(Nr,Nr)*(r/5+4/5)-A1^2*exp(1i*(-2*kx*X)) - A2^2*exp(1i*(2*kx*X)) + A3^2*exp(1i*(-2*ky*Y)) + A4^2*exp(1i*(2*ky*Y));
%     end
%     P0=E0;
% else
%     %%%%%% начальное распределение - однородный профиль
% 
    E0=ones(1,Nr)*sqrt(r_max-1-(delta./(1+sigma_min)).^2);%+0.001*randn(Nr,Nr);
    P0=E0*(1-1i*delta./(1+sigma_min));
    D0=ones(1,Nr)*(1+(delta./(1+sigma_min)).^2);
    E0(1,1)=E0(1,1)+1e-6;
% end

Ecur=E0;
Pcur=P0;
Dcur=D0;
%%
ID=int16(0);
mkdir('time');
mkdir('time/1d');
picsDir='pics';
mkdir(picsDir);
filename=[apertureType,' delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist(['time/1d/INT ',filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=[apertureType,' delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
erowfilename=['time/1d/INT ',filename,'.dat'];
fid = fopen(erowfilename, 'w');
fclose(fid);
irowfilename=['time/1d/IE ',filename,'.dat'];
fid = fopen(irowfilename, 'w');
fclose(fid);
datafilenameK=['time/1d/k ',filename,'.dat'];
wsfilename=[filename,'.mat'];
%%% custom palette
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';
            
%%
tic
% Dot=zeros(Nr,Nr);
if t==0, 
    erow=zeros(nl,1); 
    irow=zeros(nl,1); 
end;
    exp1=exp(-gamma*tau/2);
    exp2=exp(-(1+1i*delta)*tau/2);
    exp3=exp(-(sigma+1i*a*K)*tau/2);
    exp4=sigma./(1i*a*K+sigma-1-1i*delta);
    jj=0;
    F=[];
while round(t/tau)<1e8%6e4
    
    Dfl=exp1.*Dcur;
    
    Pfl=exp2.*Pcur;
    
    Ef=fft(Ecur);Pf=fft2(Pcur).*exp4;
    Efl=ifft(exp2.*Pf+exp3.*(Ef-Pf));

    %%%%%%%%%%%%%%%%% START Trapezional scheme
    E=Ecur;   P=Pcur;  D=Dcur;
    Ej=Ecur;
    exp5=Dcur.*Ecur;
    exp6=2*r-0.5*(conj(Ecur).*Pcur+Ecur.*conj(Pcur));
    Ef=fft(Efl);
    nmax=100;eps=1e-15; % Iteration Constants
	%%%% Start of Iteration Process
    for j=1:nmax
        %%%% Nonlinear part dt
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
	%%%% End of Iteration Process
    if j==jj
    else jj=j;disp(['j=',num2str(j)]);
    end;
    if j==nmax
        disp(['Not converges to ',num2str(eps),' in ',num2str(nmax),' iterations']);%pause;
    end
    Ecur=E; Pcur=P; Dcur=D;
    
    nstep=nstep+int32(1);
    t=t+tau;
    if abs(mod(nstep,int32(erow_incr)))<tau/2,
%         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
%         disp(t);
        nrec=nrec+1;
        erow(nrec)=Ecur(Nr/2);%Ecur(round(Nr/2),round(Nr/2))];
        irow(nrec)=sum(sum(abs(Ecur)));
        [m,mInd]=max(Ef);
        qmax(nrec)=mInd;
        if exist(stopfilename,'file')==2,
                break;
        end
        if nrec==nl, 
    %         Erow=[Erow erow];erow=[];%DI=[DI dI];dI=[]; % Imax=[Imax imax];imax=[];
    %         disp('write center');
            disp(['t= ',num2str(t),' at ',num2str(toc),' s']);
            fid = fopen(erowfilename, 'a');
            fwrite(fid, abs(erow).^2, 'double');
            fclose(fid); 
            fid = fopen(irowfilename, 'a');
            fwrite(fid, irow, 'double');
            fclose(fid); 
            fid = fopen(datafilenameK, 'a');
            fwrite(fid, qmax, 'double');
            fclose(fid); 
            nrec=int32(0);   % запись матрицы в файл (40 байт)
            
            %%% near field images
%             Ecur=ifft2(Ef);
%             imagesc(x0(Nr*NrecSize+1:Nr*(1-NrecSize)),x0(Nr*NrecSize+1:Nr*(1-NrecSize)),abs(Ecur(Nr*NrecSize+1:Nr*(1-NrecSize),Nr*NrecSize+1:Nr*(1-NrecSize))).^2);
%             colorbar('vert');colormap(gray);
%             xlabel('x','FontSize',14,'FontWeight','bold');ylabel('y','FontSize',14,'FontWeight','bold');
%             gcffilename=[filename,' nf ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
%             colormap(grayCustom);
%             gcffilename=[filename,' nf inv ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
%             %%% far field images
%             imagesc(abs(fftshift(Ef)));colorbar('vert');colormap(gray);
%             xlabel('k_x','FontSize',14,'FontWeight','bold');ylabel('k_y','FontSize',14,'FontWeight','bold');
%             gcffilename=[filename,' ff ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
%             colormap(grayCustom);
%             gcffilename=[filename,' ff inv ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
            %%% closing gcf
            close gcf;
            if mod(round(t/tau),1e5)<tau/2
                save(wsfilename);
            end
            KK(:,nk)=fftshift(abs(Ef).^2);
            nk=nk+1;
        end
        if Ecur(Nr/2) ~= Ecur(Nr/2)
            disp('Nan solution');
            break;
%         else
%             if sum(sum(abs(e1))) < 1e-20
%                 disp('zero solution');
%                 break;
%             end
        end
    end
%     imagesc(abs(e1).^2);colorbar('vert');
%     imagesc(abs(fftshift(Ef)).^2);colorbar('vert');
%     title(t);
%     drawnow;
end
toc;
save(wsfilename);
%%
if nrec>0, 
%     Erow=[Erow erow];erow=[];  DI=[DI dI];dI=[]; Imax=[Imax imax];imax=[]; 
%     disp('write end');
    fid = fopen(erowfilename, 'a');
    fwrite(fid, abs(erow(1:nrec)).^2, 'double');
    fclose(fid); 
    fid = fopen(irowfilename, 'a');
    fwrite(fid, irow(1:nrec), 'double');
    fclose(fid);
    fid = fopen(datafilenameK, 'a');
    fwrite(fid, qmax, 'double');
    fclose(fid); 
%     disp(nrec);
    nrec=int32(0);
    erow=zeros(nl,1);
    irow=zeros(nl,1);
    qmax=zeros(nl,1);
end;

%% LOCAL INTENSITY
fid = fopen(erowfilename, 'rb');
INT = fread(fid, 'double'); 
fclose(fid); 
% Ist=r_max-1-(delta/(1+sigma_min))^2;
figure;plot((1:numel(INT))*erow_incr*tau,INT,'black');
% title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
gcffilename=[picsDir,'/Int(t) 1D ',filename,'.jpg'];
xlabel('t','FontSize',18,'FontWeight','bold');ylabel('I=|E|^2','FontSize',18,'FontWeight','bold');
saveas(gcf,gcffilename,'jpg');%close gcf;

% INTEGRAL INTENSITY
fid = fopen(irowfilename, 'rb');
IE = fread(fid, 'double'); 
fclose(fid); 
figure;plot((1:numel(IE))*erow_incr*tau,IE.^2/Nr^4,'black');
% title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
gcffilename=[picsDir,'/Integral Int(t) 1D ',filename,'.jpg'];
xlabel('t','FontSize',18,'FontWeight','bold');ylabel('\int|E|^2','FontSize',18,'FontWeight','bold');
saveas(gcf,gcffilename,'jpg');close gcf;
%% NEAR FIELD
x0=(0:Nr-1)*h;
Ecur=ifft2(Ef);
figure;
% imagesc(x0(Nr*NrecSize+1:Nr*(1-NrecSize)),x0(Nr*NrecSize+1:Nr*(1-NrecSize)),abs(Ecur(Nr*NrecSize+1:Nr*(1-NrecSize),Nr*NrecSize+1:Nr*(1-NrecSize))).^2);
plot(x0,abs(Ecur).^2,'black','LineWidth',1);
xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
gcffilename=[picsDir,'/nearfield 1D ',filename,'.jpg'];
% colormap(grayCustom);
saveas(gcf,gcffilename,'jpg');%close gcf;
plot(x0,angle(Ecur),'black','LineWidth',1);
xlabel('x','FontSize',14,'FontWeight','bold');ylabel('y','FontSize',14,'FontWeight','bold');
gcffilename=[picsDir,'/phase 1D ',filename,'.jpg'];
saveas(gcf,gcffilename,'jpg');%close gcf;

%% FAR FIELD (STANDART)
figure;
plot(fftshift(kx),abs(fftshift(Ef)),'black','LineWidth',1);

%% FAR FIELD (CUSTOM)
% x0=(0:Nr-1)*h;
xdz=0.03;   
% x=xdz*x0+2;
x=xdz*fftshift(kx);    
Fur=abs(ft(Ecur.')); 
figure;   
% imagesc(2*pi*x(Nr*NrecSize+1:Nr*(1-NrecSize)),2*pi*x(Nr*NrecSize+1:Nr*(1-NrecSize)),Fur(Nr*3/16:Nr*13/16,Nr*3/16:Nr*13/16)); 
plot(2*pi*x,Fur,'black','LineWidth',1); 
% mesh(Fur(Nr*3/16:Nr*13/16,Nr*3/16:Nr*13/16));
xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
dzfilename=[picsDir,'/dz 1d I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% saveas(gcf,dzfilename,'png');%close gcf;

%%
% fur1=Fur;
% Kcut=7;
% for i=1:Nr
%     for j=1:Nr
%         if sqrt((2*pi*x(i))^2+(2*pi*x(j))^2)<Kcut
%             fur1(i,j)=0;
%         end
%     end
% end
% imagesc(2*pi*x,2*pi*x,(fur1));
% xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
% colorbar('vert'); colormap(gray);
% dzfilename=['dz cut',num2str(Kcut),' I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% % saveas(gcf,dzfilename,'png');%close gcf;
%%
zoomOffset=floor(t/1e2)*1e2;
zoomInc1=5e2;
zoomInc2=1e2;
figure;plot((1:numel(IE))*erow_incr*tau,IE.^2/Nr^4,'black');
xlabel('t','FontSize',18,'FontWeight','bold');
ylabel('\int|E|^2','FontSize',18,'FontWeight','bold');
xlim([zoomOffset-zoomInc1 zoomOffset])
gcffilename=['Integral Int(t) 2D ',filename,'.jpg'];
saveas(gcf,[picsDir,'/zoom ',gcffilename],'jpg');
xlim([zoomOffset-zoomInc2 zoomOffset])
saveas(gcf,[picsDir,'/zoom1 ',gcffilename],'jpg');close gcf;

figure;plot((1:numel(INT))*erow_incr*tau,INT,'black');
xlabel('t','FontSize',18,'FontWeight','bold');
ylabel('I=|E|^2','FontSize',18,'FontWeight','bold');
xlim([zoomOffset-zoomInc1 zoomOffset])
gcffilename=['Int(t) 2D ',filename,'.jpg'];
saveas(gcf,[picsDir,'/zoom ',gcffilename],'jpg');
xlim([zoomOffset-zoomInc2 zoomOffset])
saveas(gcf,[picsDir,'/zoom1 ',gcffilename],'jpg');close gcf;

%%
nw=8192;W=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);
plot(W,log(fftshift(abs(fft(INT(numel(INT)-nw+1:numel(INT)))))),'black')
xlim([0 min([10,max(W)])]);
xlabel('\omega','FontSize',18,'FontWeight','bold');
ylabel('lg(C)','FontSize',18,'FontWeight','bold');
saveas(gcf,[picsDir,'/spectra  ',filename,'.jpg'],'jpg');close gcf;

imagesc((1:size(KK,2))*tau*1e4,fftshift(kx),KK)
ylim([-25 25])
ylabel('kx','FontSize',18,'FontWeight','bold');
xlabel('t','FontSize',18,'FontWeight','bold');
colormap(grayCustom);
filename=[picsDir,'/K(t) 2D delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.jpg'];
saveas(gcf,filename,'jpg');

fid = fopen(datafilenameK, 'rb');
qindexes = fread(fid, 'double'); 
fclose(fid); 
qind=qindexes(1:10:numel(qindexes));
plot((1:numel(qind))*tau*10,kx(qind+1),'black')
xlabel('t','FontSize',18,'FontWeight','bold');
ylabel('(k_x)_m_a_x','FontSize',18,'FontWeight','bold');
ylim([-25 25]);
filename=['k(t) zoom 2D delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.jpg'];
saveas(gcf,[picsDir,'/',filename],'jpg');
xlim([numel(qind)*0.995*tau*10 numel(qind)*tau*10]);
saveas(gcf,[picsDir,'/zoom ',filename],'jpg');