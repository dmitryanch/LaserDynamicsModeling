%%%% Simple script for MB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are can be set as on infinite or on finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON delta_aperture_ON sigma_aperture_ON
T=4e5;
t=0;
tau=0.5;
Nt=ceil(T/tau);
sigma=0.025;
gamma=5e-5;
delta=-0.2;
beta=0.1;
a=6e-4;
r=2;

%%%% indicates parameter profiling
delta_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;

%%%%% Transverse shape of active medium
apertureType = 'square';
% apertureType = 'circle';

diag=1;    %% use diagonal simmetry of square vortex latticces or 
Npolos=0;
L=40;
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
h=L/128;
Nr=ceil(L/h);
stopfilename='stop0';        %%% file name that will stop calculating
kfil=int32((Nr-1)/7*3)+1;     %%% parameter for dealising

kx=((-Nr/2):(Nr/2-1))*2*pi/L;
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
nl=1000;


r_max=r;sigma_min=sigma;

% Npolos=8;disp(['?????????? ????? ?? ????????: ',num2str(Npolos)]);
% a=delta*(2*pi*Npolos/L)^-2;

d=L;
dx=h;
% x0=linspace(-d/2,d/2,Nr);
x0=(0:Nr-1)*h;
[X Y]=meshgrid(x0);
Kt=sqrt(delta/a);Kwave=0;
%%
nw=8192;w=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);
disp(['Kx_max=',num2str(max(kx))]);
disp(['dk=',num2str(kx(2))]);
disp(['W_max=',num2str(max(w))]);
disp(['dW=',num2str(2*pi/(tau*nw*erow_incr))]);
disp(['OMEGA_0_teor=',num2str(sigma*delta/(1+sigma))]);

%%
E0=0.3*rand(Nr,Nr);
P0=0.4*rand(Nr,Nr);
D0=0.5*rand(Nr,Nr);
% E0=smooth(0.001*randn(Nr,Nr),1);
% P0=smooth(0.001*randn(Nr,Nr),1);
% D0=smooth(0.001*randn(Nr,Nr),1);

% if delta >= 0
%     %%%%%% ????????? ????????????? - ?????????? ??????? ??????
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
%     %%%%%% ????????? ????????????? - ?????????? ???????
% 
%     E0=sqrt(r_max-1-(delta./(1+sigma_min)).^2).*ones(Nr,Nr);%+0.001*randn(Nr,Nr);
%     P0=E0*(1-1i*delta./(1+sigma_min));
%     D0=(1+(delta./(1+sigma_min)).^2).*ones(Nr,Nr);
%     E0=E0+1e-6*randn(Nr,Nr);
% %     E0=E0+1e-8;
% end

Pf=fft2(P0);
Ef=fft2(E0);
Df=fft2(D0);
Ef(1,1)=Ef(1,1)+1e-4;
Df=makeReal2d(Df);
Ef=filter2d(Ef);
Pf=filter2d(Pf);
Df=filter2d(Df);
%%
ID=int16(0);
mkdir('time');
mkdir('time/2d');
picsDir='pics/2d';
mkdir(picsDir);
filename=[apertureType,' delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist(['time/2d/INT ',filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=[apertureType,' delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
erowfilename=['time/2d/INT ',filename,'.dat'];
fid = fopen(erowfilename, 'w');
fclose(fid);
irowfilename=['time/2d/IE ',filename,'.dat'];
fid = fopen(irowfilename, 'w');
fclose(fid);
rerowfilename=['time/2d/Re ',filename,'.dat'];
fid = fopen(rerowfilename, 'w');
fclose(fid);
wsfilename=['2d ',filename,'.mat'];
%%% custom palette
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';
            
%%
tic
% Dot=zeros(Nr,Nr);
if t==0, 
    erow=zeros(nl,1); 
    irow=zeros(nl,1); 
    rerow=zeros(nl,1); 
end;

%%%%% precomp
    
[expM11,expM12,expM22,expM33,exp2M11,exp2M12,exp2M22,exp2M33,invM11,invM12,invM22,invM33,a211,a212,a222,a233,b211,b212,b222,b233,c111,c112,c122,c133,c211,c212,c222,c233,c311,c312,c322,c333]=precompSmbInf(sigma,gamma,delta,a,K,tau);

%%%%%
while round(t/tau)<10e4/tau%6e4
    
    [FEf,FPf,FDf,e1]=nonlinear2dInfFeedback(Ef,Pf,Df,delta,gamma,r,sigma,beta);
    
    AEf=exp2M11.*Ef+exp2M12.*Pf+a211.*FEf+a212.*FPf;
    APf=exp2M22.*Pf+a222.*FPf;
    ADf=exp2M33.*Df+a233.*FDf;
    
    [FAEf,FAPf,FADf]=nonlinear2dInfFeedback(AEf,APf,ADf,delta,gamma,r,sigma,beta);
    
    BEf=expM11.*Ef+expM12.*Pf+b211.*(2*FAEf-FEf)+b212.*(2*FAPf-FPf);
    BPf=expM22.*Pf+b222.*(2*FAPf-FPf);
    BDf=expM33.*Df+b233.*(2*FADf-FDf);
    
    [FBEf,FBPf,FBDf]=nonlinear2dInfFeedback(BEf,BPf,BDf,delta,gamma,r,sigma,beta);
    
    cEf=expM11.*Ef+expM12.*Pf+c111.*FEf+c112.*FPf+4*c211.*FAEf+4*c212.*FAPf+c311.*FBEf+c312.*FBPf;
    cPf=expM22.*Pf+c122.*FPf+4*c222.*FAPf+c322.*FBPf;
    cDf=expM33.*Df+c133.*FDf+4*c233.*FADf+c333.*FBDf;
    
    
    Ef=filter2d(cEf);
    Pf=filter2d(cPf);
    Df=filter2d(cDf);
    
    nstep=nstep+int32(1);
    t=t+tau;
    if abs(mod(nstep,int32(erow_incr)))<tau/2,
%         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
%         disp(t);
        nrec=nrec+1;
        erow(nrec)=e1(Nr/2,Nr/2);%Ecur(round(Nr/2),round(Nr/2))];
        irow(nrec)=sum(sum(abs(e1)));
        rerow(nrec)=real(e1(Nr/2,Nr/2));
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
            fid = fopen(rerowfilename, 'a');
            fwrite(fid, rerow, 'double');
            fclose(fid);
            nrec=int32(0);   % ?????? ??????? ? ???? (40 ????)
            
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
%             close gcf;
%             save(wsfilename);
        end
        if e1(Nr/2,Nr/2) ~= e1(Nr/2,Nr/2)
            disp('Nan solution');
            break;
%         else
%             if sum(sum(abs(e1))) < 1e-20
%                 disp('zero solution');
%                 break;
%             end
        end
    end
    imagesc(abs(e1).^2);colorbar('vert');
%     imagesc(abs(fftshift(Ef)).^2);colorbar('vert');
%     title(t);
    drawnow;
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
    fid = fopen(rerowfilename, 'a');
    fwrite(fid, rerow(1:nrec), 'double');
    fclose(fid);
%     disp(nrec);
    nrec=int32(0);
    erow=zeros(nl,1);
    irow=zeros(nl,1);
    rerow=zeros(nl,1);
end;

%% LOCAL INTENSITY
fid = fopen(erowfilename, 'rb');
INT = fread(fid, 'double'); 
fclose(fid); 
% Ist=r_max-1-(delta/(1+sigma_min))^2;
figure;plot((1:numel(INT))*erow_incr*tau,INT,'black');
% title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
gcffilename=[picsDir,'/Int(t) 2D ',filename,'.jpg'];
xlabel('t','FontSize',18,'FontWeight','bold');ylabel('I=|E|^2','FontSize',18,'FontWeight','bold');
% saveas(gcf,gcffilename,'jpg');%close gcf;

% INTEGRAL INTENSITY
fid = fopen(irowfilename, 'rb');
IE = fread(fid, 'double'); 
fclose(fid); 
figure;plot((1:numel(IE))*erow_incr*tau,IE.^2/Nr^4,'black');
% title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
gcffilename=[picsDir,'/Integral Int(t) 2D ',filename,'.jpg'];
xlabel('t','FontSize',18,'FontWeight','bold');ylabel('\int|E|^2','FontSize',18,'FontWeight','bold');
% saveas(gcf,gcffilename,'jpg');close gcf;
%%
% REAL PART OF FIELD
fid = fopen(rerowfilename, 'rb');
RE = fread(fid, 'double'); 
fclose(fid); 
figure;plot((1:numel(RE))*erow_incr*tau,RE,'black');
% title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
gcffilename=[picsDir,'/re(E(t)) 1D ',filename,'.jpg'];
xlabel('t','FontSize',18,'FontWeight','bold');ylabel('Re(E)','FontSize',18,'FontWeight','bold');
% saveas(gcf,gcffilename,'jpg');%close gcf;
%%
nw=min(8192/4,numel(RE));%
plot((numel(RE)-nw+1:numel(RE))*erow_incr*tau,RE(numel(RE)-nw+1:numel(RE)),'black');
W=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);
figure;plot(W,((log(fftshift(abs(fft(RE(numel(RE)-nw+1:numel(RE)))).^2)))),'black')
xlim([0 min([10,max(W)])]);
xlabel('\omega','FontSize',18,'FontWeight','bold');
ylabel('lg(ReE)','FontSize',18,'FontWeight','bold');
% saveas(gcf,[picsDir,'/spectra RE3  ',filename,'.jpg'],'jpg');%close gcf;
%% NEAR FIELD
x0=(0:Nr-1)*h;
Ecur=ifft2(Ef);
figure;
% imagesc(x0(Nr*NrecSize+1:Nr*(1-NrecSize)),x0(Nr*NrecSize+1:Nr*(1-NrecSize)),abs(Ecur(Nr*NrecSize+1:Nr*(1-NrecSize),Nr*NrecSize+1:Nr*(1-NrecSize))).^2);
imagesc(x0,x0,abs(Ecur).^2);
xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
gcffilename=[picsDir,'/nearfield 2D ',filename,'.jpg'];
% colormap(grayCustom);
colormap(gray);
colorbar('vert');
% saveas(gcf,gcffilename,'jpg');%close gcf;
% figure;imagesc(x0,x0,angle(Ecur));colorbar('vert');
xlabel('x','FontSize',14,'FontWeight','bold');ylabel('y','FontSize',14,'FontWeight','bold');
gcffilename=[picsDir,'/phase 2D ',filename,'.jpg'];
% saveas(gcf,gcffilename,'jpg');%close gcf;

%% FAR FIELD (STANDART)
figure;
dzmult=1;
imagesc(dzmult*fftshift(kx),dzmult*fftshift(kx),(sqrt(abs(fftshift(Ef)))));
lim=30;
xlim([-lim lim]);ylim([-lim lim]);
xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
colorbar('vert'); colormap(gray);
dzfilename=[picsDir,'/dz1 I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% saveas(gcf,dzfilename,'png');%close gcf;
%%
absE=(abs(Ef).^2).^0.5;
sumE=sum(sum(absE));
[Yy, KK, Ynorm]=getDistribByWaveNumber2D(absE,kx,128);
hold on;grid on;plot(KK,Yy/sumE*100,'k','LineWidth',2);hold off;
% hold on;grid on;plot(KK,Ynorm/sum(Ynorm)*100,'k','LineWidth',2);hold off;
xlim([0 min([15 max(KK)])]);
set(gca,'FontSize',fontSize);xlabel('k','FontSize',labelSize,'FontWeight','bold');ylabel('%','FontSize',labelSize,'FontWeight','bold');
saveas(gcf,['NrgByK ',filename,' ',num2str(t),'.png'],'png');

%%
ff=abs(fftshift(Ef)).^2;
allField = sum(sum(ff));
ff1=ff;
Kcut=12;
kx1=fftshift(kx);
for i=1:Nr
    for j=1:Nr
        if sqrt(kx1(i)^2+kx1(j)^2)<=Kcut
            ff1(i,j)=0;
        end
    end
end
cutField=sum(sum(ff1));
imagesc(fftshift(kx),fftshift(kx),ff1);
title(['cut=',num2str(Kcut),' => ',num2str((allField-cutField)/allField*100),'%'])
%% FAR FIELD (CUSTOM)
% x0=(0:Nr-1)*h;
xdz=0.05;   
% x=xdz*x0+2;
x=xdz*fftshift(kx);    
Ecur=ifft2(Ef);
Fur=dz(Ecur); 
figure;   
% imagesc(2*pi*x(Nr*NrecSize+1:Nr*(1-NrecSize)),2*pi*x(Nr*NrecSize+1:Nr*(1-NrecSize)),Fur(Nr*3/16:Nr*13/16,Nr*3/16:Nr*13/16)); 
imagesc(2*pi*x,2*pi*x,Fur); 
% mesh(Fur(Nr*3/16:Nr*13/16,Nr*3/16:Nr*13/16));
xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
colorbar('vert'); colormap(gray);
dzfilename=[picsDir,'/dz I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% saveas(gcf,dzfilename,'png');%close gcf;

%%
fur1=Fur;
Kcut=5;
for i=1:Nr
    for j=1:Nr
        if sqrt((2*pi*x(i))^2+(2*pi*x(j))^2)<Kcut
            fur1(i,j)=0;
        end
    end
end
imagesc(2*pi*x,2*pi*x,(fur1));
xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
colorbar('vert'); colormap(gray);
dzfilename=[picsDir,'/dz cut',num2str(Kcut),' I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
saveas(gcf,dzfilename,'png');%close gcf;
%%
x0=(0:Nr-1)*h;
Ecur=ifft2(Ef);
mesh(x0,x0,abs(Ecur).^2);colormap(gray);
zlim([0 1.2*max(max(abs(Ecur).^2))])
xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',14,'FontWeight','bold');     
zlabel('I=|E|^2','FontSize',18,'FontWeight','bold')
saveas(gcf,[picsDir,'/mesh ',filename,'.png'],'png');%close gcf;
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
nw=min(8192/4,numel(INT));W=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);
plot((numel(INT)-nw+1:numel(INT))*erow_incr*tau,INT(numel(INT)-nw+1:numel(INT)),'black');
figure;plot(W,log((fftshift(abs(fft(INT(numel(INT)-nw+1:numel(INT)))))).^0.1),'black')
xlim([0 min([1.5 max(W)])]);
xlabel('\omega','FontSize',18,'FontWeight','bold');
ylabel('lg(C)','FontSize',18,'FontWeight','bold');
% saveas(gcf,[picsDir,'/spectra  ',filename,'.jpg'],'jpg');%close gcf;