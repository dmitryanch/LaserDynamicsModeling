%%%% Simple script for VMB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are tuned as infinite and finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON sigma_aperture_ON a K
T=4e5;
t=0;
tau=5e-3;
Nt=ceil(T/tau);

r=5;
omega=1;
a=.01;
sigma=0.1;
b=0.01;
c=0.05;

Ne = 7; % number of equations
omega_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;
lyapunov_ON = false;
k0=real(sqrt(omega/a));
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
% Ep0=smooth(n0*randn(Nr,Nr),1);
% Em0=smooth(n0*randn(Nr,Nr),1);
% Pp0=smooth(n0*randn(Nr,Nr),1);
% Pm0=smooth(n0*randn(Nr,Nr),1);
% Np0=smooth(n0*randn(Nr,Nr),1);
% Nm0=smooth(n0*randn(Nr,Nr),1);
% M0 =smooth(n0*randn(Nr,Nr),1);

%%%%%% Single Travelling co-waves
Ep0=sqrt(2*b*c/(3*c+b)*(r-1))*exp(1i*(k0/sqrt(2)*(X+Y)));
Em0=sqrt(2*b*c/(3*c+b)*(r-1))*exp(1i*(k0/sqrt(2)*(X+Y)));
Pp0=Ep0;
Pm0=Em0;
Np0=abs(Ep0).^2*3/2/b;
Nm0=abs(Ep0).^2*3/2/b;
M0=abs(Ep0).^2/2/c;

%%%%%%  Two Travelling trans-waves
% Ep0=sqrt(2*b*c/(3*c+b)*(r-1))*(exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y))));
% Em0=sqrt(2*b*c/(3*c+b)*(r-1))*(exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y))));
% Pp0=Ep0;
% Pm0=Em0;
% Np0=abs(Ep0).^2*3/b*(3+exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y)))+0.5*exp(1i*(k0/sqrt(2)*(X-Y)))+0.5*exp(-1i*(k0/sqrt(2)*(X-Y))));
% Nm0=abs(Ep0).^2*3/b*(3+0.5*exp(1i*(k0/sqrt(2)*(X+Y)))+0.5*exp(-1i*(k0/sqrt(2)*(X+Y)))+exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y))));
% M0=abs(Ep0).^2/2/c*(exp(1i*(k0/sqrt(2)*(X)))+exp(-1i*(k0/sqrt(2)*(Y)))+exp(1i*(k0/sqrt(2)*(-X)))+exp(-1i*(k0/sqrt(2)*(-Y))));

%%%%%%  4SVL
Ep00=sqrt((2*b*c*r*sigma)/(5*(gammaA + sigma)*(b + 3*c)) - (2*b*c*(a*k0^2 + gammaP - omega)^2)/(5*(b + 3*c)*(gammaA + sigma + 1)^2) - (2*b*c)/(5*(b + 3*c)));
Ep0=Ep00*(exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y)))+exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y)+pi)));
Em0=Ep0;
Pp0=Ep0*(a*k0^2*1i + gammaA + gammaP*1i + sigma - w*1i)/sigma;
Pm0=Em0*(a*k0^2*1i + gammaA + gammaP*1i + sigma - w*1i)/sigma;
Np0=Ep00^2*(6/b+6*gammaA/(b*sigma)+exp(-k0/sqrt(2)*(X+Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma))-exp(k0/sqrt(2)*(X+Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma))+exp(-k0/sqrt(2)*(X-Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma))-exp(k0/sqrt(2)*(X-Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma)));
Nm0=Np0;
M0=Ep00^2*(2/c+2*gammaA/(c*sigma)+exp(-k0/sqrt(2)*(X+Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma))-exp(k0/sqrt(2)*(X+Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma))+exp(-k0/sqrt(2)*(X-Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma))-exp(k0/sqrt(2)*(X-Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma)));

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
nl=1000;

SP=0;
SPar=[];
    
ID=int16(0);
filename=['INT2D omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' gA=',num2str(gammaA),' gP=',num2str(gammaP),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist([filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=['INT2D omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' gA=',num2str(gammaA),' gP=',num2str(gammaP),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
datafilename=[filename,'.dat'];
while exist(['INT2Dp ',filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' gA=',num2str(gammaA),' gP=',num2str(gammaP),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
end
datafilenameP=['INT2Dp ',filename,'.dat'];
fid = fopen(datafilenameP, 'w');
fclose(fid);
datafilenameM=['INT2Dm ',filename,'.dat'];
fid = fopen(datafilenameM, 'w');
fclose(fid);
energyfilenameP = ['Energy2DP ',filename,'.dat'];
energyfilenameM = ['Energy2DM ',filename,'.dat'];
%%% custom palette
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';
% init=zeros(Ne,Nr,Nr);
%     init(1,:,:)=Epf;
%     init(2,:,:)=Emf;
%     init(3,:,:)=Ppf;
%     init(4,:,:)=Pmf;
%     init(5,:,:)=Npf;
%     init(6,:,:)=Nmf;
%     init(7,:,:)=Mf;
    
% A=zeros(Ne,Nr,Nr);
% B=zeros(Ne,Nr,Nr);
% C=zeros(Ne,Nr,Nr);
DifF=[];
Dif=[];


%%%%% matrix's precomp
[expM11 expM13 expM22 expM24 expM31 expM33 expM42 expM44 expM55 expM66 expM77 exp2M11 exp2M13 exp2M22 exp2M24 exp2M31 exp2M33 exp2M42 exp2M44 exp2M55 exp2M66 exp2M77 invM11 invM13 invM22 invM24 invM31 invM33 invM42 invM44 invM55 invM66 invM77 a211 a213 a222 a224 a231 a233 a242 a244 a255 a266 a277 b211 b213 b222 b224 b231 b233 b242 b244 b255 b266 b277 c111 c113 c122 c124 c131 c133 c142 c144 c155 c166 c177 c211 c213 c222 c224 c231 c233 c242 c244 c255 c266 c277 c311 c313 c322 c324 c331 c333 c342 c344 c355 c366 c377] = precomp(K,tau,a,omega,sigma,r,b,c);
%%%%%
            
%%
ffttime=0;
nltime=0;
tstart=tic
% Dot=zeros(Nr,Nr);
if t==0, 
    erowP=zeros(nl,1); 
    erowM=zeros(nl,1); 
    energyP=zeros(nl,1);
    energyM=zeros(nl,1); 
end;

while round(t/tau)*tau<1e7*tau%6e4
    
    %%% matrix form
%     init(1,:,:)=Epf;
%     init(2,:,:)=Emf;
%     init(3,:,:)=Ppf;
%     init(4,:,:)=Pmf;
%     init(5,:,:)=Npf;
%     init(6,:,:)=Nmf;
%     init(7,:,:)=Mf;
%     ntime=tic;
    [fEpf,fEmf,fPpf,fPmf,fNpf,fNmf,fMf, ep, em]=nonlin2d(Epf,Emf,Ppf,Pmf,Npf,Nmf,Mf);
%     nltime=nltime+toc(ntime);
%     ffttime=ffttime+f1;
    aEpf=exp2M11.*Epf+exp2M13.*Ppf+a211.*fEpf+a213.*fPpf;
    aEmf=exp2M22.*Emf+exp2M24.*Pmf+a222.*fEmf+a224.*fPmf;
    aPpf=exp2M31.*Epf+exp2M33.*Ppf+a231.*fEpf+a233.*fPpf;
    aPmf=exp2M42.*Emf+exp2M44.*Pmf+a242.*fEmf+a244.*fPmf;
	aNpf=exp2M55.*Npf+a255.*fNpf;
	aNmf=exp2M66.*Nmf+a266.*fNmf;
	aMf=exp2M77.*Mf+a277.*fMf;
    
%     ntime=tic;
    [faEpf,faEmf,faPpf,faPmf,faNpf,faNmf,faMf]=nonlin2d(aEpf,aEmf,aPpf,aPmf,aNpf,aNmf,aMf);
%     nltime=nltime+toc(ntime);
%     ffttime=ffttime+f1;
	bEpf=expM11.*Epf+expM13.*Ppf+b211.*(2*faEpf-fEpf)+b213.*(2*faPpf-fPpf);
    bEmf=expM22.*Emf+expM24.*Pmf+b222.*(2*faEmf-fEmf)+b224.*(2*faPmf-fPmf);
    bPpf=expM31.*Epf+expM33.*Ppf+b231.*(2*faEpf-fEpf)+b233.*(2*faPpf-fPpf);
    bPmf=expM42.*Emf+expM44.*Pmf+b242.*(2*faEmf-fEmf)+b244.*(2*faPmf-fPmf);
	bNpf=expM55.*Npf+b255.*(2*faNpf-fNpf);
	bNmf=expM66.*Nmf+b266.*(2*faNmf-fNmf);
	bMf=expM77.*Mf+b277.*(2*faMf-fMf);
%     ntime=tic;
    [fbEpf,fbEmf,fbPpf,fbPmf,fbNpf,fbNmf,fbMf]=nonlin2d(bEpf,bEmf,bPpf,bPmf,bNpf,bNmf,bMf);
%     nltime=nltime+toc(ntime);
%     ffttime=ffttime+f1;
    cEpf=expM11.*Epf+expM13.*Ppf+c111.*fEpf+c113.*fPpf+4*c211.*faEpf+4*c213.*faPpf+c311.*fbEpf+c313.*fbPpf;
    cEmf=expM22.*Emf+expM24.*Pmf+c122.*fEmf+c124.*fPmf+4*c222.*faEmf+4*c224.*faPmf+c322.*fbEmf+c324.*fbPmf;
    cPpf=expM31.*Epf+expM33.*Ppf+c131.*fEpf+c133.*fPpf+4*c231.*faEpf+4*c233.*faPpf+c331.*fbEpf+c333.*fbPpf;
    cPmf=expM42.*Emf+expM44.*Pmf+c142.*fEmf+c144.*fPmf+4*c242.*faEmf+4*c244.*faPmf+c342.*fbEmf+c344.*fbPmf;
	cNpf=expM55.*Npf+c155.*fNpf+4*c255.*faNpf+c355.*fbNpf;
	cNmf=expM66.*Nmf+c166.*fNmf+4*c266.*faNmf+c366.*fbNmf;
	cMf=expM77.*Mf+c177.*fMf+4*c277.*faMf+c377.*fbMf;
    
    Epf=filter2d(cEpf);
    Emf=filter2d(cEmf);
    Ppf=filter2d(cPpf);
    Pmf=filter2d(cPmf);
    Npf=filter2d(cNpf);
    Nmf=filter2d(cNmf);
    Mf=filter2d(cMf);
    
    nstep=nstep+int32(1);
    t=t+tau;
    if abs(mod(nstep,int32(erow_incr)))<tau/2,
%         dI=[dI (max(max(abs(Ecur).^2))-min(min(abs(Ecur).^2)))/max(max(abs(Ecur).^2))];
%         disp(t);

        if(ep(1) ~= ep(1))
            disp('NAN detected');
            break;
        end
        nrec=nrec+1;
        erowP(nrec)=abs(ep(Nr/2,Nr/2))^2;
        erowM(nrec)=abs(em(Nr/2,Nr/2))^2;%Ecur(round(Nr/2),round(Nr/2))];
        energyP(nrec)=sum(sum(abs(ep).^2));
        energyM(nrec)=sum(sum(abs(em).^2));
        
        if nrec==nl, 
    %         Erow=[Erow erow];erow=[];%DI=[DI dI];dI=[]; % Imax=[Imax imax];imax=[];
    %         disp('write center');
            disp(['t= ',num2str(t),' at ',num2str(toc(tstart)),' s']);
            disp(['sumI= ',num2str((energyP(nrec)+energyM(nrec))/Nr^2)]);
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
            
            figure;colormap(grayCustom);
            imagesc(abs(ep));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,['prof/ep ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(em));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,['prof/em ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(ep) + abs(em));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,['prof/e ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(fftshift(Epf)));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,['prof/epf ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(fftshift(Emf)));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,['prof/emf ',filename,' ',num2str(t),'.png'],'png');
            imagesc(fftshift(abs(Epf) + abs(Emf)));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,['prof/ef ',filename,' ',num2str(t),'.png'],'png');
            close gcf;
            save(['2D inf ',filename,'.mat']);
            
            if exist(stopfilename,'file')==2,
                break;
            end
        end
%         imagesc(abs(ep).^2+abs(em).^2);colorbar('vert');
%         drawnow;
    end
    
% % %     plot(1:Nr,abs(ep).^2+abs(em).^2);
% % %     title(t);
% % %     drawnow;
end
toc(tstart);

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
plot((1:numel(ENERGYp))*erow_incr*tau,ENERGYp./Nr^2);
fid = fopen(energyfilenameM, 'rb');
ENERGYm = fread(fid, 'double'); 
fclose(fid); 
subplot('224');
plot((1:numel(ENERGYm))*erow_incr*tau,ENERGYm./Nr^2);
gcffilename=['2Dtotal omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.jpg'];
xlabel('Время t');ylabel('Уровень интенсивности I=|E|^2');%saveas(gcf,gcffilename,'jpg');close gcf;

%%
save(['2D ',filename,'.mat']);