%%%% Simple script for VMB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are tuned as infinite and finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON sigma_aperture_ON a K
T=4e5;
t=0;
tau=1e-2;
Nt=ceil(T/tau);

r=10;
omega=2;
a=.01;
sigma=0.01;
b=0.0001;
c=0.005;

Ne = 7; % number of equations
omega_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;
lyapunov_ON = false;

diag=1;
Npolos=20;
k0=real(sqrt(omega/a));
if omega > 0
    L=2*pi/k0*Npolos;
else 
    L = 1;
end

h=L/256;
Nr=ceil(L/h);
stopfilename='stop';
kfil=int32((Nr-1)/7*3)+1;

r_max=r;sigma_min=sigma;
thr=1+omega^2/(1+sigma)^2*(1-sign(omega))/2;
r=r_max*thr;

l=2;
d=L;
dx=h;
x0=(0:Nr-1)*h;
[X Y]=meshgrid(x0);

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
% n0=0.001;
% Ep0=(n0*randn(1,Nr));
% Em0=(n0*randn(1,Nr));
% Pp0=(n0*randn(1,Nr));
% Pm0=(n0*randn(1,Nr));
% Np0=(n0*randn(1,Nr));
% Nm0=(n0*randn(1,Nr));
% M0 =(n0*randn(1,Nr));

%%%%%% Single Travelling co-waves
Ep0=sqrt(2*b*c/(3*c+b)*(r-1))*exp(1i*(k0*x0));
Em0=sqrt(2*b*c/(3*c+b)*(r-1))*exp(1i*(k0*x0));
Pp0=Ep0;
Pm0=Em0;
Np0=abs(Ep0).^2*3/2/b;
Nm0=abs(Ep0).^2*3/2/b;
M0=abs(Ep0).^2/2/c;

%%%%%% Single Travelling counter-waves 
% Ep0=sqrt(2*b*c/(3*c+b)*(r-1))*exp(1i*(k0*x0));
% Em0=sqrt(2*b*c/(3*c+b)*(r-1))*exp(-1i*(k0*x0));
% Pp0=Ep0;
% Pm0=Em0;
% Np0=abs(Ep0).^2*3/2/b;
% Nm0=abs(Ep0).^2*3/2/b;
% M0=abs(Ep0).^2/2/c.*exp(2*1i*k0*x0);

%%%%%%% Standing-waves solution
% I0=(2*b*c*(r - 1))/(3*b + 9*c);
% Ep0=sqrt(I0)*exp(1i*(k0*x0))+sqrt(I0)*exp(-1i*(k0*x0));
% Em0=sqrt(I0)*exp(1i*(k0*x0))+sqrt(I0)*exp(-1i*(k0*x0));
% Pp0=Ep0;
% Pm0=Em0;
% Np0=(3*exp(-k0*x0*2i).*I0.*(exp(k0*x0*2i) + 1).^2)/(2*b);
% Nm0=Np0;
% M0=(exp(-k0*x0*2i)*I0.*(exp(k0*x0*2i) + 1).^2)/(2*c);

Epf=fft(Ep0);
Emf=fft(Em0);
Ppf=fft(Pp0);
Pmf=fft(Pm0);
Npf=fft(Np0);
Nmf=fft(Nm0);
Mf =fft(M0);
Epf(1,1)=Epf(1,1)+1e-4;
Npf=makeReal1d(Npf);
Nmf=makeReal1d(Nmf);
Epf=filter1d(Epf);
Emf=filter1d(Emf);
Ppf=filter1d(Ppf);
Pmf=filter1d(Pmf);
Npf=filter1d(Npf);
Nmf=filter1d(Nmf);
Mf=filter1d(Mf);

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
filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist(['INTp ',filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
datafilenameP=['INTp ',filename,'.dat'];
fid = fopen(datafilenameP, 'w');
fclose(fid);
datafilenameM=['INTm ',filename,'.dat'];
fid = fopen(datafilenameM, 'w');
fclose(fid);
energyfilenameP = ['EnergyP ',filename,'.dat'];
energyfilenameM = ['EnergyM ',filename,'.dat'];
%%% custom palette
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';

%%
%%%%% matrix's precomp
cc=zeros(Ne,Ne,Nr);
ecc=zeros(Ne,Ne,Nr);
ecc2=zeros(Ne,Ne,Nr);
aecc=zeros(Ne,Ne,Nr);
becc=zeros(Ne,Ne,Nr);
for q=1:Nr
    cc(:,:,q)=[-sigma-1i*a*K(q),0,sigma,0,0,0,0;
        0,-sigma-1i*a*K(q),0,sigma,0,0,0;
        r,0,-1-1i*omega,0,0,0,0;
        0,r,0,-1-1i*omega,0,0,0;
        0,0,0,0,-b,0,0;
        0,0,0,0,0,-b,0;
        0,0,0,0,0,0,-c];

    ecc(:,:,q)=expm(cc(:,:,q)*tau);
    ecc2(:,:,q)=expm(cc(:,:,q)*0.5*tau);
    
    aecc(:,:,q)=(ecc2(:,:,q)-eye(Ne))/cc(:,:,q);
    becc(:,:,q)=(ecc(:,:,q)-eye(Ne))/cc(:,:,q);
end
%     invc=inv(c);
ttau=tau;
ccc1=zeros(Ne,Ne,Nr);ccc2=zeros(Ne,Ne,Nr);ccc3=zeros(Ne,Ne,Nr);
for q=1:Nr
    ccc1(:,:,q) = -4/ttau^2*(cc(:,:,q)^-1)^3 - 1/ttau*(cc(:,:,q)^-1)^2 + ecc(:,:,q)*4/ttau^2*(cc(:,:,q)^-1)^3 - ecc(:,:,q)*3/ttau*(cc(:,:,q)^-1)^2 + ecc(:,:,q)*(cc(:,:,q)^-1);
    ccc2(:,:,q) = 2/ttau^2*(cc(:,:,q)^-1)^3 + 1/ttau*(cc(:,:,q)^-1)^2 + ecc(:,:,q)*(-2)/ttau^2*(cc(:,:,q)^-1)^3 + ecc(:,:,q)/ttau*(cc(:,:,q)^-1)^2;
    ccc3(:,:,q) = -4/ttau^2*(cc(:,:,q)^-1)^3 - 3/ttau*(cc(:,:,q)^-1)^2 - cc(:,:,q)^-1 + ecc(:,:,q)*4/ttau^2*(cc(:,:,q)^-1)^3 - ecc(:,:,q)/ttau*(cc(:,:,q)^-1)^2;
    
%     ccc1(:,:,q)=(-4*eye(Ne)-ttau*cc(:,:,q)+ecc(:,:,q)*(4*eye(Ne)-3*ttau*cc(:,:,q)+ttau^2*cc(:,:,q)^2))/ttau^2/(cc(:,:,q)^3);
%     ccc2(:,:,q)=(2*eye(Ne)+ttau*cc(:,:,q)+ecc(:,:,q)*(-2*eye(Ne)+ttau*cc(:,:,q)))/ttau^2/(cc(:,:,q)^3);
%     ccc3(:,:,q)=(-4*eye(Ne)-3*ttau*cc(:,:,q)-ttau^2*cc(:,:,q)^2+ecc(:,:,q)*(4*eye(Ne)-ttau*cc(:,:,q)))/ttau^2/(cc(:,:,q)^3);
end
    
%%%%%
A=zeros(Ne,Nr);
B=zeros(Ne,Nr);
C=zeros(Ne,Nr);
DifF=[];
Dif=[];
%%
tic
if t==0, 
    erowP=zeros(nl,1); 
    erowM=zeros(nl,1); 
    energyP=zeros(nl,1);
    energyM=zeros(nl,1);
end;
while round(t/tau)*tau<100e5*tau%6e4
    
    %%% matrix form
    init=[Epf;Emf;Ppf;Pmf;Npf;Nmf;Mf];
    [F, ep, em]=nonlinear1d(init);
    
    for q=1:Nr
        A(:,q)=ecc2(:,:,q)*init(:,q)+aecc(:,:,q)*F(:,q);
    end

    [FA]=nonlinear1d(A);
    for q=1:Nr
        B(:,q)=ecc(:,:,q)*init(:,q)+becc(:,:,q)*(2*FA(:,q)-F(:,q));
    end

    [FB]=nonlinear1d(B);
    for q=1:Nr
        C(:,q)=ecc(:,:,q)*init(:,q)+ccc1(:,:,q)*F(:,q)+4*ccc2(:,:,q)*FA(:,q)+ccc3(:,:,q)*FB(:,q);
    end

    Epf=filter1d(C(1,:));
    Emf=filter1d(C(2,:));
    Ppf=filter1d(C(3,:));
    Pmf=filter1d(C(4,:));
    Npf=filter1d(C(5,:));
    Nmf=filter1d(C(6,:));
    Mf =filter1d(C(7,:));
    
    nstep=nstep+int32(1);
    t=t+tau;
    if abs(mod(nstep,int32(erow_incr)))<tau/2,
%         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
%         disp(t);

        if(ep(1) ~= ep(1))
            break;
        end
        nrec=nrec+1;
        erowP(nrec)=abs(ep(Nr/2))^2;
        erowM(nrec)=abs(em(Nr/2))^2;%Ecur(round(Nr/2),round(Nr/2))];
        energyP(nrec)=sum(abs(ep).^2);
        energyM(nrec)=sum(abs(em).^2);
        
        
         if nrec==nl, 
    %         Erow=[Erow erow];erow=[];%DI=[DI dI];dI=[]; % Imax=[Imax imax];imax=[];
    %         disp('write center');
            disp(['t= ',num2str(t),' at ',num2str(toc),' s']);
            disp(['I= ',num2str(erowP(nrec)+erowM(nrec))]);
            fid = fopen(datafilenameP, 'a');
            fwrite(fid, erowP, 'double');
            fclose(fid); 
            fid = fopen(datafilenameM, 'a');
            fwrite(fid, erowM, 'double');
            fclose(fid); 
            
            fid = fopen(energyfilenameP, 'a');
            fwrite(fid, energyP, 'double');
            fclose(fid); 
            fid = fopen(energyfilenameM, 'a');
            fwrite(fid, energyM, 'double');
            fclose(fid); 
            nrec=int32(0);   % запись матрицы в файл (40 байт)
            if exist(stopfilename,'file')==2,
                break;
            end
            
        end
    end
    
%     plot(1:Nr,abs(ep).^2+abs(em).^2);
% %     subplot(211);    
% %     plot(1:Nr,abs(ep).^2);
% %     subplot(212);
% %     plot(1:Nr,abs(em).^2);
% %     title(t);
% %     drawnow;
end
toc;

%%
if nrec>0, 
%     Erow=[Erow erow];erow=[];  DI=[DI dI];dI=[]; Imax=[Imax imax];imax=[]; 
%     disp('write end');
    fid = fopen(datafilenameP, 'a');
    fwrite(fid, erowP(1:nrec), 'double');
    fclose(fid); 
    fid = fopen(datafilenameM, 'a');
    fwrite(fid, erowM(1:nrec), 'double');
    fclose(fid); 
    
    fid = fopen(energyfilenameP, 'a');
    fwrite(fid, energyP(1:nrec), 'double');
    fclose(fid); 
    fid = fopen(energyfilenameM, 'a');
    fwrite(fid, energyM(1:nrec), 'double');
    fclose(fid); 
%     disp(nrec);
    nrec=int32(0);
    erowP=zeros(nl,1);
    erowM=zeros(nl,1);
    energyP=zeros(nl,1);
    energyM=zeros(nl,1);
end;


%%
fid = fopen(datafilenameP, 'rb');
INTp = fread(fid, 'double'); 
fclose(fid); 
subplot('221');
plot((1:numel(INTp))*erow_incr*tau,INTp);
title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a)]);
fid = fopen(datafilenameM, 'rb');
INTm = fread(fid, 'double'); 
fclose(fid); 
subplot('222');
plot((1:numel(INTm))*erow_incr*tau,INTm);
fid = fopen(energyfilenameP, 'rb');
ENERGYp = fread(fid, 'double'); 
fclose(fid); 
subplot('223');
plot((1:numel(ENERGYp))*erow_incr*tau,ENERGYp./Nr);
fid = fopen(energyfilenameM, 'rb');
ENERGYm = fread(fid, 'double'); 
fclose(fid); 
subplot('224');
plot((1:numel(ENERGYm))*erow_incr*tau,ENERGYm./Nr);


filename=['Int(t) 1D omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.jpg'];
% xlabel('Время t');ylabel('Уровень интенсивности I=|E|^2');saveas(gcf,filename,'jpg');close gcf;
