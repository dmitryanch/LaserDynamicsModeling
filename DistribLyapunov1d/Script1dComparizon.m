%%%% Simple script for MB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are tuned as infinite and finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON delta_aperture_ON sigma_aperture_ON a K
T=4e5;
t=0;
tau=1e-2;
Nt=ceil(T/tau);
sigma=1;
gamma=2.2;
delta=0;
a=.01;
r=60;

delta_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;
lyapunov_ON = true;

diag=1;
Npolos=30;
L=1;
if delta > 0
    if diag == 1
        L=2*pi/sqrt(delta/a/2)*Npolos;
    else 
        L=2*pi/sqrt(delta/a)*Npolos;
    end
else
%     L=1;
    L=2*pi/35.5884443058738*12;
end
h=L/128;
Nr=ceil(L/h);
stopfilename='stop1';
kfil=int32((Nr-1)/5)+1;
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

l=2;
d=L;
dx=h;
% x0=linspace(-d/2,d/2,Nr);
x0=(0:Nr-1)*h;
[X Y]=meshgrid(x0);
Kt=sqrt(delta/a);Kwave=0;

%%% накачка - квадрат
if r_aperture_ON==1
    r=-1*ones(1,Nr);
    r(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=r_max;
    r=smooth(r,1);  
end
if sigma_aperture_ON==1
    sigma=sigma_min*5*ones(1,Nr);
    sigma(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=sigma_min;
    sigma=smooth(sigma,1);
end

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
E0=1/sqrt(2)*ones(1,Nr).*sqrt(r_max-1-(delta./(1+sigma_min)).^2).*r/r_max+1i/sqrt(2)*ones(1,Nr).*sqrt(r_max-1-(delta./(1+sigma_min)).^2).*r/r_max;
P0=E0*(1-1i*delta./(1+sigma_min));
D0=ones(1,Nr).*(1+(delta./(1+sigma_min)).^2).*r/r_max;
E0(1,1)=E0(1,1)+1e-4;

Pf=fft(P0);
Ef=fft(E0);
Df=fft(D0);
Ef(1,1)=Ef(1,1)+1e-4;
Df=makeReal1d(Df);
Ef=filter1d(Ef);
Pf=filter1d(Pf);
Df=filter1d(Df);
ef=Ef;
pf=Pf;
df=Df;

% kx=linspace(-Nr/2,Nr/2,Nr)*2*pi/L;
kx=((-Nr/2):(Nr/2-1))*2*pi/L;
kx=fftshift(kx);
K=zeros(1,Nr);
for i=1:Nr
    K(i)=kx(i)^2;
end
nstep=int32(0);
nrec=int32(0);
erow_incr=10;
nl=1000;

SP=0;
SPar=[];
    
ID=int16(0);
filename=['INT delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist([filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=['INT delta=',num2str(delta),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' gamma=',num2str(gamma),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
datafilename=[filename,'.dat'];
%%% custom palette
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';
            
%%
tic
% Dot=zeros(Nr,Nr);
if t==0, 
    erow=zeros(nl,1); 
end;

%%%%% precomp
    
    if sigma_aperture_ON==1
        ccc = -a*1i*K-1;
    else
        ccc = -a*1i*K-sigma;
    end
    cc1(1,:)=(-4-ccc*tau+exp(ccc*tau).*(4-3*ccc*tau+(ccc*tau).^2))./(ccc.^3*tau^2);
    cc1(2,:)=(2+ccc*tau+exp(ccc*tau).*(-2+ccc*tau))./(ccc.^3*tau^2);
    cc1(3,:)=(-4-3*ccc*tau-(ccc*tau).^2+exp(ccc*tau).*(4-ccc*tau))./(ccc.^3*tau^2);

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

%%%%% Lyapunov's precomp
cc=zeros(3,3,Nr);ecc=zeros(3,3,Nr);
for q=1:Nr
    cc(:,:,q)=[-sigma-1i*a*K(q),0,0;
        0,-1-1i*delta,0;
        0,0,-gamma];
    if sigma_aperture_ON==1
        ccclr = -1;
        cc(1,1,q)=-1i*a*K(q)-1;
    end
    ecc(:,:,q)=expm(cc(:,:,q)*tau);
    ecc2(:,:,q)=expm(cc(:,:,q)*0.5*tau);
end
%     invc=inv(c);
    ttau=tau;
    ccc1=zeros(3,3,Nr);ccc2=zeros(3,3,Nr);ccc3=zeros(3,3,Nr);
    for q=1:Nr
        ccc1(:,:,q)=(-4*eye(3)-ttau*cc(:,:,q)+ecc(:,:,q)*(4*eye(3)-3*ttau*cc(:,:,q)+ttau^2*cc(:,:,q)^2))/ttau^2/(cc(:,:,q)^3);
        ccc2(:,:,q)=(2*eye(3)+ttau*cc(:,:,q)+ecc(:,:,q)*(-2*eye(3)+ttau*cc(:,:,q)))/ttau^2/(cc(:,:,q)^3);
        ccc3(:,:,q)=(-4*eye(3)-3*ttau*cc(:,:,q)-ttau^2*cc(:,:,q)^2+ecc(:,:,q)*(4*eye(3)-ttau*cc(:,:,q)))/ttau^2/(cc(:,:,q)^3);
    end
    
%%%%%
A=zeros(3,Nr);
B=zeros(3,Nr);
C=zeros(3,Nr);
DifF=[];
Dif=[];
while round(t/tau)*tau<100e5*tau%6e4
    
    [FEf,FPf,FDf,e1]=nonlinear1d(Ef,Pf,Df,sigma,delta,gamma,r);
    
    AEf=Ef.*exp(0.5*ccc*tau)+(exp(0.5*ccc*tau)-1).*FEf./ccc;
    APf=Pf.*exp(0.5*ddd*tau)+(exp(0.5*ddd*tau)-1).*FPf./ddd;
    ADf=Df.*exp(0.5*eee*tau)+(exp(0.5*eee*tau)-1).*FDf./eee;
    
    [FAEf,FAPf,FADf]=nonlinear1d(AEf,APf,ADf,sigma,delta,gamma,r);
    
    BEf=Ef.*exp(ccc*tau)+(exp(ccc*tau)-1).*(2*FAEf-FEf)./ccc;
    BPf=Pf.*exp(ddd*tau)+(exp(ddd*tau)-1).*(2*FAPf-FPf)./ddd;
    BDf=Df.*exp(eee*tau)+(exp(eee*tau)-1).*(2*FADf-FDf)./eee;
    
    [FBEf,FBPf,FBDf]=nonlinear1d(BEf,BPf,BDf,sigma,delta,gamma,r);
    
    Ef=Ef.*exp(ccc*tau)+cc1(1,:).*FEf +4*FAEf.*cc1(2,:)+FBEf.*cc1(3,:);
    Pf=Pf.*exp(ddd*tau)+cc2(1).*FPf     +4*FAPf.*cc2(2)    +FBPf.*cc2(3);
    Df=Df.*exp(eee*tau)+cc3(1).*FDf     +4*FADf.*cc3(2)    +FBDf.*cc3(3);
    
    Ef=filter1d(Ef);
    Pf=filter1d(Pf);
    Df=filter1d(Df);
    
    %%% matrix form
    init=[ef;pf;df];
    [Fef,Fpf,Fdf,e11]=nonlinear1d(ef,pf,df,sigma,delta,gamma,r);
    F=[Fef; Fpf; Fdf];

    for q=1:Nr
        A(:,q)=ecc2(:,:,q)*init(:,q)+(ecc2(:,:,q)-eye(3))/cc(:,:,q)*F(:,q);
    end

    [FAef, FApf, FAdf]=nonlinear1d(A(1,:), A(2,:), A(3,:),sigma,delta,gamma,r);
    FA=[FAef; FApf; FAdf];
    for q=1:Nr
        B(:,q)=ecc(:,:,q)*init(:,q)+(ecc(:,:,q)-eye(3))/cc(:,:,q)*(2*FA(:,q)-F(:,q));
    end

    [FBef, FBpf, FBdf]=nonlinear1d(B(1,:), B(2,:), B(3,:),sigma,delta,gamma,r);
    FB=[FBef; FBpf; FBdf];
    for q=1:Nr
        C(:,q)=ecc(:,:,q)*init(:,q)+ccc1(:,:,q)*F(:,q)+4*ccc2(:,:,q)*FA(:,q)+ccc3(:,:,q)*FB(:,q);
    end

    ef=filter1d(C(1,:));
    pf=filter1d(C(2,:));
    df=filter1d(C(3,:));
    
    nstep=nstep+int32(1);
    t=t+tau;
    if abs(mod(nstep,int32(erow_incr)))<tau/2,
%         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
%         disp(t);
        nrec=nrec+1;
        erow(nrec)=e1(Nr/2);%Ecur(round(Nr/2),round(Nr/2))];
        
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
    
        
          
        if abs(mod(nstep,int32(20)))<tau/2,
            difF=max(abs(Ef-ef));
            dif=max(abs(ifft(Ef)-ifft(ef)));
            disp(['t= ',num2str(t),' difF= ',num2str(difF),' dif=',num2str(dif)]);
            DifF=[DifF difF];
            Dif=[Dif dif];
            subplot(211);
            plot(linspace(0,t,numel(DifF)),DifF);

            subplot(212);
            plot(linspace(0,t,numel(Dif)),Dif);
            drawnow;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     imagesc(abs(e1).^2);colorbar('vert');
%     imagesc(abs(fftshift(Ef)).^2);colorbar('vert');
%     title(t);
%     drawnow;
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
% xlabel('Время t');ylabel('Уровень интенсивности I=|E|^2');saveas(gcf,filename,'jpg');close gcf;
