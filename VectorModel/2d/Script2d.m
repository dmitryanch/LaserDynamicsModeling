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
sigma=0.1;
b=0.01;
c=0.05;

Ne = 7; % number of equations
omega_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;
lyapunov_ON = false;

diag=1;
Npolos=4;
if omega > 0
    L=2*pi/sqrt(omega/a/2)*Npolos;
else 
    L = 1;
end

h=L/128;
Nr=ceil(L/h);
stopfilename='stop1';
kfil=int32((Nr-1)/7*3)+1;

r_max=r;sigma_min=sigma;

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
n0=0.1;
Ep0=smooth(n0*randn(Nr,Nr),1);
Em0=smooth(n0*randn(Nr,Nr),1);
Pp0=smooth(n0*randn(Nr,Nr),1);
Pm0=smooth(n0*randn(Nr,Nr),1);
Np0=smooth(n0*randn(Nr,Nr),1);
Nm0=smooth(n0*randn(Nr,Nr),1);
M0 =smooth(n0*randn(Nr,Nr),1);

%%%%%% начальное распределение - однородный профиль
% E0=1/sqrt(2)*ones(1,Nr).*sqrt(r_max-1-(delta./(1+sigma_min)).^2).*r/r_max+1i/sqrt(2)*ones(1,Nr).*sqrt(r_max-1-(delta./(1+sigma_min)).^2).*r/r_max;
% P0=E0*(1-1i*delta./(1+sigma_min));
% D0=ones(1,Nr).*(1+(delta./(1+sigma_min)).^2).*r/r_max;
% E0(1,1)=E0(1,1)+1e-4;

Epf=fft2(Ep0);
Emf=fft2(Em0);
Ppf=fft2(Pp0);
Pmf=fft2(Pm0);
Npf=fft2(Np0);
Nmf=fft2(Nm0);
Mf =fft2(M0);
Epf(1,1)=Epf(1,1)+1e-4;
Npf=makeReal2d(Npf);
Nmf=makeReal2d(Nmf);
Epf=filter2d(Epf);
Emf=filter2d(Emf);
Ppf=filter2d(Ppf);
Pmf=filter2d(Pmf);
Npf=filter2d(Npf);
Nmf=filter2d(Nmf);
Mf=filter2d(Mf);

% kx=linspace(-Nr/2,Nr/2,Nr)*2*pi/L;
kx=((-Nr/2):(Nr/2-1))*2*pi/L;
kx=fftshift(kx);
K=zeros(Nr,Nr);
for i=1:Nr
    for j=1:Nr
        K(i,j)=kx(i)^2+kx(j)^2;
    end
end
nstep=int32(0);
nrec=int32(0);
erow_incr=10;
nl=10;

SP=0;
SPar=[];
    
ID=int16(0);
filename=['INT2D omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist([filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=['INT2D omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
datafilename=[filename,'.dat'];
%%% custom palette
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';
init=zeros(Ne,Nr,Nr);
    init(1,:,:)=Epf;
    init(2,:,:)=Emf;
    init(3,:,:)=Ppf;
    init(4,:,:)=Pmf;
    init(5,:,:)=Npf;
    init(6,:,:)=Nmf;
    init(7,:,:)=Mf;
    
A=zeros(Ne,Nr,Nr);
B=zeros(Ne,Nr,Nr);
C=zeros(Ne,Nr,Nr);
DifF=[];
Dif=[];


%%%%% matrix's precomp
cc=zeros(Ne,Ne,Nr,Nr);ecc=zeros(Ne,Ne,Nr,Nr);ecc2=zeros(Ne,Ne,Nr,Nr);
for p=1:Nr
    for q=1:Nr
        cc(:,:,p,q)=[-sigma-1i*a*K(p,q),0,sigma,0,0,0,0;
            0,-sigma-1i*a*K(p,q),0,sigma,0,0,0;
            r,0,-1-1i*omega,0,0,0,0;
            0,r,0,-1-1i*omega,0,0,0;
            0,0,0,0,-b,0,0;
            0,0,0,0,0,-b,0;
            0,0,0,0,0,0,-c];

        ecc(:,:,p,q)=expm(cc(:,:,p,q)*tau);
        ecc2(:,:,p,q)=expm(cc(:,:,p,q)*0.5*tau);
    end
end
%     invc=inv(c);
ttau=tau;
ccc1=zeros(Ne,Ne,Nr,Nr);ccc2=zeros(Ne,Ne,Nr,Nr);ccc3=zeros(Ne,Ne,Nr,Nr);
for p=1:Nr
    for q=1:Nr
%         ccc1(:,:,p,q)=(-4*eye(Ne)-ttau*cc(:,:,p,q)+ecc(:,:,p,q)*(4*eye(Ne)-3*ttau*cc(:,:,p,q)+ttau^2*cc(:,:,p,q)^2))/ttau^2/cc(:,:,p,q)/cc(:,:,p,q)/cc(:,:,p,q);
%         ccc2(:,:,p,q)=(2*eye(Ne)+ttau*cc(:,:,p,q)+ecc(:,:,p,q)*(-2*eye(Ne)+ttau*cc(:,:,p,q)))/ttau^2/cc(:,:,p,q)/cc(:,:,p,q)/cc(:,:,p,q);
%         ccc3(:,:,p,q)=(-4*eye(Ne)-3*ttau*cc(:,:,p,q)-ttau^2*cc(:,:,p,q)^2+ecc(:,:,p,q)*(4*eye(Ne)-ttau*cc(:,:,p,q)))/ttau^2/cc(:,:,p,q)/cc(:,:,p,q)/cc(:,:,p,q);
        
        ccc1(:,:,p,q) = -4/ttau^2*(cc(:,:,p,q)^-1)^3 - 1/ttau*(cc(:,:,p,q)^-1)^2 + ecc(:,:,p,q)*4/ttau^2*(cc(:,:,p,q)^-1)^3 - ecc(:,:,p,q)*3/ttau*(cc(:,:,p,q)^-1)^2 + ecc(:,:,p,q)*(cc(:,:,p,q)^-1);
        ccc2(:,:,p,q) = 2/ttau^2*(cc(:,:,p,q)^-1)^3 + 1/ttau*(cc(:,:,p,q)^-1)^2 + ecc(:,:,p,q)*(-2)/ttau^2*(cc(:,:,p,q)^-1)^3 + ecc(:,:,p,q)/ttau*(cc(:,:,p,q)^-1)^2;
        ccc3(:,:,p,q) = -4/ttau^2*(cc(:,:,p,q)^-1)^3 - 3/ttau*(cc(:,:,p,q)^-1)^2 - cc(:,:,p,q)^-1 + ecc(:,:,p,q)*4/ttau^2*(cc(:,:,p,q)^-1)^3 - ecc(:,:,p,q)/ttau*(cc(:,:,p,q)^-1)^2;
    
    end
end
    
%%%%%
            
%%
ffttime=0;
nltime=0;
tstart=tic
% Dot=zeros(Nr,Nr);
if t==0, 
    erow=zeros(nl,1); 
end;

while round(t/tau)*tau<1e2*tau%6e4
    
    %%% matrix form
%     init(1,:,:)=Epf;
%     init(2,:,:)=Emf;
%     init(3,:,:)=Ppf;
%     init(4,:,:)=Pmf;
%     init(5,:,:)=Npf;
%     init(6,:,:)=Nmf;
%     init(7,:,:)=Mf;
    ntime=tic;
    [F, ep, em,f1]=nonlinear2d(init);
    nltime=nltime+toc(ntime);
    ffttime=ffttime+f1;
    for p=1:Nr
        for q=1:Nr
            A(:,p,q)=ecc2(:,:,p,q)*init(:,p,q)+(ecc2(:,:,p,q)-eye(Ne))/cc(:,:,p,q)*F(:,p,q);
        end
    end
    ntime=tic;
    [FA,ep,em,f1]=nonlinear2d(A);
    nltime=nltime+toc(ntime);
    ffttime=ffttime+f1;
    for p=1:Nr
        for q=1:Nr
            B(:,p,q)=ecc(:,:,p,q)*init(:,p,q)+(ecc(:,:,p,q)-eye(Ne))/cc(:,:,p,q)*(2*FA(:,p,q)-F(:,p,q));
        end
    end
    ntime=tic;
    [FB,ep,em,f1]=nonlinear2d(B);
    nltime=nltime+toc(ntime);
    ffttime=ffttime+f1;
    for p=1:Nr
        for q=1:Nr
            C(:,p,q)=ecc(:,:,p,q)*init(:,p,q)+ccc1(:,:,p,q)*F(:,p,q)+4*ccc2(:,:,p,q)*FA(:,p,q)+ccc3(:,:,p,q)*FB(:,p,q);
        end
    end
    
    

    init(1,:,:)=filter2d(C(1,:,:));
    init(2,:,:)=filter2d(C(2,:,:));
    init(3,:,:)=filter2d(C(3,:,:));
    init(4,:,:)=filter2d(C(4,:,:));
    init(5,:,:)=filter2d(C(5,:,:));
    init(6,:,:)=filter2d(C(6,:,:));
    init(7,:,:)=filter2d(C(7,:,:));
    
    nstep=nstep+int32(1);
    t=t+tau;
%     if abs(mod(nstep,int32(erow_incr)))<tau/2,
% %         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
% %         disp(t);
% 
%         if(ep(1) ~= ep(1))
%             break;
%         end
%         nrec=nrec+1;
%         erow(nrec)=abs(ep(Nr/2))^2+abs(em(Nr/2))^2;%Ecur(round(Nr/2),round(Nr/2))];
%         
%         if nrec==nl, 
%     %         Erow=[Erow erow];erow=[];%DI=[DI dI];dI=[]; % Imax=[Imax imax];imax=[];
%     %         disp('write center');
%             disp(['t= ',num2str(t),' at ',num2str(toc),' s']);
%             disp(['I= ',num2str(erow(nrec))]);
%             fid = fopen(datafilename, 'a');
%             fwrite(fid, erow, 'double');
%             fclose(fid); 
%             nrec=int32(0);   % запись матрицы в файл (40 байт)
%             if exist(stopfilename,'file')==2,
%                 break;
%             end
%         end
%         imagesc(abs(ep).^2+abs(em).^2);
%         drawnow;
%     end
    
%     plot(1:Nr,abs(ep).^2+abs(em).^2);
%     title(t);
%     drawnow;
end
toc(tstart);

%%
if nrec>0, 
%     Erow=[Erow erow];erow=[];  DI=[DI dI];dI=[]; Imax=[Imax imax];imax=[]; 
%     disp('write end');
    fid = fopen(datafilename, 'a');
    fwrite(fid, erow(1:nrec), 'double');
    fclose(fid); 
%     disp(nrec);
    nrec=int32(0);
    erow=zeros(nl,1);
end;


%%
% INT=[];
% fid = fopen(datafilename, 'rb');
% INT = fread(fid, 'double'); 
% fclose(fid); 
% figure;plot((1:numel(INT))*erow_incr*tau,INT);
% title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a)]);
% filename=['Int(t) 2D omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.jpg'];
% xlabel('Время t');ylabel('Уровень интенсивности I=|E|^2');saveas(gcf,filename,'jpg');close gcf;
