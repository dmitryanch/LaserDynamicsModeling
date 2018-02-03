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
sigma=0.1;
gamma=0.001;
delta=-0.25;
a=.01;
r=5;
m=0.6;
wRel=sqrt(2*gamma*sigma*(r-1-delta^2)/(1+delta^2));
w=2*wRel;

%%%% indicates parameter profiling
delta_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;

%%%%% Transverse shape of active medium
apertureType = 'square';
% apertureType = 'circle';

diag=1;    %% use diagonal simmetry of square vortex latticces or 
Npolos=0;
L=20;
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
h=L/256;
Nr=ceil(L/h);
stopfilename='stop0';        %%% file name that will stop calculating
kfil=int32((Nr-1)/7*3)+1;     %%% parameter for dealising
% int Nr = 256, Nt=1;
%     double t = 0.0, tau = 0.01, L = 5.0, h = L/Nr;		
%     double tau_2 = tau / 2;
%     double gamma =0.1;
%     double sigma =10;
%     double delta = -4.5;
%     double a =0.01;
%     double r = 100.0;

r_max=r;sigma_min=sigma;

% Npolos=8;disp(['?????????? ????? ?? ????????: ',num2str(Npolos)]);
% a=delta*(2*pi*Npolos/L)^-2;

d=L;
dx=h;
% x0=linspace(-d/2,d/2,Nr);
x0=(0:Nr-1)*h;
[X Y]=meshgrid(x0);
Kt=sqrt(delta/a);Kwave=0;

r=ones(Nr,Nr)*r_max;
if r_aperture_ON==1
    r=-1*ones(Nr,Nr);
    if apertureType == 'square'
        %%% ??????? - ???????
        r(round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)),round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)))=r_max;
    end
    if apertureType == 'circle'
        %%%% ??????? - ????
        for m=1:Nr 
            for n=1:Nr
                if (m-Nr/2-0.5)^2+(n-Nr/2-0.5)^2<=Nr^2*l^2/4
                    r(m,n)=r_max;  
                end
            end
        end
    end

end

sigma=sigma_min*ones(Nr,Nr);       
if sigma_aperture_ON==1
    sigma=sigma_min*1e2*ones(Nr,Nr);        
    if apertureType == 'square'
    %%% ??????? - ???????
        sigma(round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)),round(Nr*0.5*(1-l)):round(Nr*0.5*(1+l)))=sigma_min;
    end
    if apertureType == 'circle'
        %%%% ??????? - ????
        for m=1:Nr 
            for n=1:Nr
                if (m-Nr/2-0.5)^2+(n-Nr/2-0.5)^2<=Nr^2*l^2/4
                    sigma(m,n)=sigma_min; 
                end
            end
        end
    end
end
r=smooth(r,1);  
sigma=smooth(sigma,1);

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
nl=1000;

nw=8192;W=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);

disp(['Kx_max=',num2str(max(kx))]);
disp(['dk=',num2str(kx(2))]);
disp(['W_max=',num2str(max(W))]);
disp(['dW=',num2str(2*pi/(tau*nw*erow_incr))]);
disp(['OMEGA_0_teor=',num2str(sigma_min*delta/(1+sigma_min))]);

%%
% E0=0.3*rand(Nr,Nr);
% P0=0.4*rand(Nr,Nr);
% D0=0.5*rand(Nr,Nr);
% % E0=smooth(0.001*randn(Nr,Nr),1);
% % P0=smooth(0.001*randn(Nr,Nr),1);
% % D0=smooth(0.001*randn(Nr,Nr),1);

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
    E0=sqrt(r_max-1-(delta./(1+sigma_min)).^2).*r/r_max;%+0.001*randn(Nr,Nr);
    P0=E0*(1-1i*delta./(1+sigma_min));
    D0=(1+(delta./(1+sigma_min)).^2).*r/r_max;
    E0(1,1)=E0(1,1)+1e-6;
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
filename=[apertureType,' delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist(['INT ',filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=[apertureType,' delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
erowfilename=['INT ',filename,'.dat'];
fid = fopen(erowfilename, 'w');
fclose(fid);
irowfilename=['IE ',filename,'.dat'];
fid = fopen(irowfilename, 'w');
fclose(fid);
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

%%%%% precomp
    
    if sigma_aperture_ON==1
        ccc = -a*1i*K-1;
    else
        ccc = -a*1i*K-sigma;
    end
    cc1(:,:,1)=(-4-ccc*tau+exp(ccc*tau).*(4-3*ccc*tau+(ccc*tau).^2))./(ccc.^3*tau^2);
    cc1(:,:,2)=(2+ccc*tau+exp(ccc*tau).*(-2+ccc*tau))./(ccc.^3*tau^2);
    cc1(:,:,3)=(-4-3*ccc*tau-(ccc*tau).^2+exp(ccc*tau).*(4-ccc*tau))./(ccc.^3*tau^2);

    ddd=-1.0;
    if delta_aperture_ON==0
        ddd = ddd - 1i*delta;
    end
    cc2(1)    =(-4-ddd*tau+exp(ddd*tau).*(4-3*ddd*tau+(ddd*tau).^2))./(ddd.^3*tau^2);
    cc2(2)    =(2+ddd*tau+exp(ddd*tau).*(-2+ddd*tau))./(ddd.^3*tau^2);
    cc2(3)    =(-4-3*ddd*tau-(ddd*tau).^2+exp(ddd*tau).*(4-ddd*tau))./(ddd.^3*tau^2);
  
    eee=-gamma;
    cc3(1)	  =(-4-eee*tau+exp(eee*tau).*(4-3*eee*tau+(eee*tau).^2))./(eee.^3*tau^2);
    cc3(2)    =(2+eee*tau+exp(eee*tau).*(-2+eee*tau))./(eee.^3*tau.^2);
    cc3(3)    =(-4-3*eee*tau-(eee*tau).^2+exp(eee*tau).*(4-eee*tau))./(eee.^3*tau^2);

%%%%%

while round(t/tau)<5e5%6e4
    r=r_max*(1+m*sin(w*t));
    [FEf,FPf,FDf,e1]=nonlinear2d(Ef,Pf,Df,sigma,delta,gamma,r);
    
    AEf=Ef.*exp(0.5*ccc*tau)+(exp(0.5*ccc*tau)-1).*FEf./ccc;
    APf=Pf.*exp(0.5*ddd*tau)+(exp(0.5*ddd*tau)-1).*FPf./ddd;
    ADf=Df.*exp(0.5*eee*tau)+(exp(0.5*eee*tau)-1).*FDf./eee;
    
    [FAEf,FAPf,FADf]=nonlinear2d(AEf,APf,ADf,sigma,delta,gamma,r);
    
    BEf=Ef.*exp(ccc*tau)+(exp(ccc*tau)-1).*(2*FAEf-FEf)./ccc;
    BPf=Pf.*exp(ddd*tau)+(exp(ddd*tau)-1).*(2*FAPf-FPf)./ddd;
    BDf=Df.*exp(eee*tau)+(exp(eee*tau)-1).*(2*FADf-FDf)./eee;
    
    [FBEf,FBPf,FBDf]=nonlinear2d(BEf,BPf,BDf,sigma,delta,gamma,r);
    
    Ef=Ef.*exp(ccc*tau)+cc1(:,:,1).*FEf +4*FAEf.*cc1(:,:,2)+FBEf.*cc1(:,:,3);
    Pf=Pf.*exp(ddd*tau)+cc2(1).*FPf     +4*FAPf.*cc2(2)    +FBPf.*cc2(3);
    Df=Df.*exp(eee*tau)+cc3(1).*FDf     +4*FADf.*cc3(2)    +FBDf.*cc3(3);
    
    Ef=filter2d(Ef);
    Pf=filter2d(Pf);
    Df=filter2d(Df);
    
    nstep=nstep+int32(1);
    t=t+tau;
    if abs(mod(nstep,int32(erow_incr)))<tau/2,
%         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
%         disp(t);
        nrec=nrec+1;
        erow(nrec)=e1(Nr/2,Nr/2);%Ecur(round(Nr/2),round(Nr/2))];
        irow(nrec)=sum(sum(abs(e1)));
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
            close gcf;
            save(wsfilename);
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
%     disp(nrec);
    nrec=int32(0);
    erow=zeros(nl,1);
    irow=zeros(nl,1);
end;

%% LOCAL INTENSITY
fid = fopen(erowfilename, 'rb');
INT = fread(fid, 'double'); 
fclose(fid); 
% Ist=r_max-1-(delta/(1+sigma_min))^2;
figure;plot((1:numel(INT))*erow_incr*tau,INT,'black');
% title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
gcffilename=['Int(t) 2D ',filename,'.jpg'];
xlabel('t','FontSize',18,'FontWeight','bold');ylabel('I=|E|^2','FontSize',18,'FontWeight','bold');
saveas(gcf,gcffilename,'jpg');%close gcf;

% INTEGRAL INTENSITY
fid = fopen(irowfilename, 'rb');
IE = fread(fid, 'double'); 
fclose(fid); 
figure;plot((1:numel(IE))*erow_incr*tau,IE.^2/Nr^4,'black');
% title(['delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a)]);
gcffilename=['Integral Int(t) 2D ',filename,'.jpg'];
xlabel('t','FontSize',18,'FontWeight','bold');ylabel('\int|E|^2','FontSize',18,'FontWeight','bold');
saveas(gcf,gcffilename,'jpg');close gcf;
%% NEAR FIELD
x0=(0:Nr-1)*h;
Ecur=ifft2(Ef);
figure;
% imagesc(x0(Nr*NrecSize+1:Nr*(1-NrecSize)),x0(Nr*NrecSize+1:Nr*(1-NrecSize)),abs(Ecur(Nr*NrecSize+1:Nr*(1-NrecSize),Nr*NrecSize+1:Nr*(1-NrecSize))).^2);
imagesc(x0,x0,abs(Ecur).^2);
xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
gcffilename=['nearfield 2D ',filename,'.jpg'];
% colormap(grayCustom);
colormap(gray);
colorbar('vert');
saveas(gcf,gcffilename,'jpg');%close gcf;
imagesc(x0,x0,angle(Ecur));colorbar('vert');
xlabel('x','FontSize',14,'FontWeight','bold');ylabel('y','FontSize',14,'FontWeight','bold');
gcffilename=['phase 2D ',filename,'.jpg'];
saveas(gcf,gcffilename,'jpg');%close gcf;

%% FAR FIELD (STANDART)
figure;
imagesc(fftshift(kx),fftshift(kx),abs(fftshift(Ef)));

%% FAR FIELD (CUSTOM)
% x0=(0:Nr-1)*h;
xdz=0.05;   
% x=xdz*x0+2;
x=xdz*fftshift(kx);    
Fur=dz(Ecur); 
figure;   
% imagesc(2*pi*x(Nr*NrecSize+1:Nr*(1-NrecSize)),2*pi*x(Nr*NrecSize+1:Nr*(1-NrecSize)),Fur(Nr*3/16:Nr*13/16,Nr*3/16:Nr*13/16)); 
imagesc(2*pi*x,2*pi*x,Fur); 
% mesh(Fur(Nr*3/16:Nr*13/16,Nr*3/16:Nr*13/16));
xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
colorbar('vert'); colormap(gray);
dzfilename=['dz I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% saveas(gcf,dzfilename,'png');%close gcf;

%%
fur1=Fur;
Kcut=0.5;
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
dzfilename=['dz cut I ',num2str(round(t/tau)),' dpi',num2str(Nr),' s',num2str(sigma_min),' r',num2str(r_max),' a',num2str(a),' d',num2str(delta),'.png'];
% saveas(gcf,dzfilename,'png');%close gcf;
%%
x0=(0:Nr-1)*h;
mesh(x0,x0,abs(Ecur).^2);colormap(gray);
xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',14,'FontWeight','bold');     
zlabel('I=|E|^2','FontSize',18,'FontWeight','bold')
saveas(gcf,['mesh ',gcffilename],'png');close gcf;

zoomOffset=7e3;
zoomInc1=1e3;
zoomInc2=1e2;
figure;plot((1:numel(IE))*erow_incr*tau,IE.^2/Nr^4,'black');
xlabel('t','FontSize',18,'FontWeight','bold');
ylabel('\int|E|^2','FontSize',18,'FontWeight','bold');
xlim([zoomOffset zoomOffset+zoomInc1])
gcffilename=['Integral Int(t) 2D ',filename,'.jpg'];
saveas(gcf,['zoom ',gcffilename],'jpg');
xlim([zoomOffset zoomOffset+zoomInc2])
saveas(gcf,['zoom1 ',gcffilename],'jpg');close gcf;

figure;plot((1:numel(INT))*erow_incr*tau,INT,'black');
xlabel('t','FontSize',18,'FontWeight','bold');
ylabel('I=|E|^2','FontSize',18,'FontWeight','bold');
xlim([zoomOffset zoomOffset+zoomInc1])
gcffilename=['Int(t) 2D ',filename,'.jpg'];
saveas(gcf,['zoom ',gcffilename],'jpg');
xlim([zoomOffset zoomOffset+zoomInc2])
saveas(gcf,['zoom1 ',gcffilename],'jpg');close gcf;

nw=8192;W=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);
plot(W,log(fftshift(abs(fft(INT(numel(INT)-nw+1:numel(INT)))))),'black')
xlim([0 10]);
xlabel('\omega','FontSize',18,'FontWeight','bold');
ylabel('lg(C)','FontSize',18,'FontWeight','bold');
saveas(gcf,['spectra  ',filename,'.jpg'],'jpg');close gcf;