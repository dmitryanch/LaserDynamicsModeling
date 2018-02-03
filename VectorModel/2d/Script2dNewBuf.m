%%%% Simple script for VMB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are tuned as infinite and finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON sigma_aperture_ON a K
T=4e5;
t=0;
tau=5e-1;
Nt=ceil(T/tau);

r=5;
omega=1;
a=.01;
sigma=0.01;
b=0.0001;
c=0.005;
gammaA=0.005;
gammaP=0.005;
dispMode(omega,a,sigma,b,c,gammaA,gammaP,r);
gcaFontSize=16;

Ne = 7; % number of equations
omega_aperture_ON = false;
sigma_aperture_ON = true;
r_aperture_ON = true;
lyapunov_ON = false;
k0=real(sqrt(omega/a));
k1=real(sqrt(1/a));
diag=1;
Npolos=17;
% if omega >= 0
    L=2*pi/k1*Npolos;
% else 
%     L = 1;
% end

h=L/256;
Nr=ceil(L/h);
stopfilename='stop1';
kfil=int32((Nr-1)/7*3)+1;

r_max=r;sigma_min=sigma;

l=1.2;
d=L;
dx=h;
x0=(0:Nr-1)*h;
[X,Y]=meshgrid(x0);
% kx=linspace(-Nr/2,Nr/2,Nr)*2*pi/L;

kx=((-Nr/2):(Nr/2-1))*2*pi/L;
kx=fftshift(kx);
K=zeros(Nr,Nr);
for i=1:Nr
    for j=1:Nr
        K(i,j)=kx(i)^2+kx(j)^2;
    end
end
%%% накачка - квадрат
if r_aperture_ON==1
    r=-1*ones(Nr,Nr);
    r(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)),round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=r_max;
    r=smooth(r,1);  
end
if sigma_aperture_ON==1
    sigma=sigma_min*1e2*ones(Nr,Nr);
    sigma(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)),round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=sigma_min;
    sigma=smooth(sigma,1);
end
% r=ones(Nr,Nr)*r_max;sigma=ones(Nr,Nr)*sigma_min;

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
% Ep0=smooth(n0*randn(Nr,Nr),1);
% Em0=smooth(n0*randn(Nr,Nr),1);
% Pp0=smooth(n0*randn(Nr,Nr),1);
% Pm0=smooth(n0*randn(Nr,Nr),1);
% Np0=smooth(n0*randn(Nr,Nr),1);
% Nm0=smooth(n0*randn(Nr,Nr),1);
% M0 =smooth(n0*randn(Nr,Nr),1);

%%%%%% Single Travelling co-waves
% Ep0=sqrt(2*b*c/(3*c+b)*(r_max-1))*exp(1i*(k0/sqrt(2)*(X+Y)));
% Em0=sqrt(2*b*c/(3*c+b)*(r_max-1))*exp(1i*(k0/sqrt(2)*(X+Y)));
% Pp0=Ep0;
% Pm0=Em0;
% Np0=abs(Ep0).^2*3/2/b;
% Nm0=abs(Ep0).^2*3/2/b;
% M0=abs(Ep0).^2/2/c;

%%%%%%  Two Travelling trans-waves
% Ep0=sqrt(2*b*c/(3*c+b)*(r-1))*(exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y))));
% Em0=sqrt(2*b*c/(3*c+b)*(r-1))*(exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y))));
% Pp0=Ep0;
% Pm0=Em0;
% Np0=abs(Ep0).^2*3/b*(3+exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y)))+0.5*exp(1i*(k0/sqrt(2)*(X-Y)))+0.5*exp(-1i*(k0/sqrt(2)*(X-Y))));
% Nm0=abs(Ep0).^2*3/b*(3+0.5*exp(1i*(k0/sqrt(2)*(X+Y)))+0.5*exp(-1i*(k0/sqrt(2)*(X+Y)))+exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y))));
% M0=abs(Ep0).^2/2/c*(exp(1i*(k0/sqrt(2)*(X)))+exp(-1i*(k0/sqrt(2)*(Y)))+exp(1i*(k0/sqrt(2)*(-X)))+exp(-1i*(k0/sqrt(2)*(-Y))));

%%%%%%  4SVL
Ep00=sqrt((2*b*c*r_max*sigma_min)/(5*(gammaA + sigma_min)*(b + 3*c)) - (2*b*c*(a*k0^2 + gammaP - omega)^2)/(5*(b + 3*c)*(gammaA + sigma_min + 1)^2) - (2*b*c)/(5*(b + 3*c)));
Ep0=Ep00*(exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y)))+exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y)+pi)));
Em0=Ep0;
Pp0=Ep0;
Pm0=Em0;
Np0=Ep00^2*(6/b+6*gammaA/(b*sigma_min)+exp(-k0/sqrt(2)*(X+Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma_min))-exp(k0/sqrt(2)*(X+Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma_min))+exp(-k0/sqrt(2)*(X-Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma_min))-exp(k0/sqrt(2)*(X-Y)*2i)*(3i/(2*b)+gammaA*3i/(2*b*sigma_min)));
Nm0=Np0;
M0=Ep00^2*(2/c+2*gammaA/(c*sigma_min)+exp(-k0/sqrt(2)*(X+Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma_min))-exp(k0/sqrt(2)*(X+Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma_min))+exp(-k0/sqrt(2)*(X-Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma_min))-exp(k0/sqrt(2)*(X-Y)*2i)*(1i/(2*c)+gammaA*1i/(2*c*sigma_min)));

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

nstep=int32(0);
nrec=int32(0);
erow_incr=10;
nl=1000;

SP=0;
SPar=[];
    
ID=int16(0);
filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP),' ID#',num2str(ID)];
while exist(['INT2Dp ',filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP),' ID#',num2str(ID)];
end
datafilename=[filename,'.dat'];
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
[expM11, expM12, expM21, expM22, expM33, expM44, expM55, expM66, expM77, exp2M11, exp2M12, exp2M21, exp2M22, exp2M33, exp2M44, exp2M55, exp2M66, exp2M77, invM11, invM12, invM21, invM22, invM33, invM44, invM55, invM66, invM77, a211, a212, a221, a222, a233, a244, a255, a266, a277, b211, b212, b221, b222, b233, b244, b255, b266, b277, c111, c112, c121, c122, c133, c144, c155, c166, c177, c211, c212, c221, c222, c233, c244, c255, c266, c277, c311, c312, c321, c322, c333, c344, c355, c366, c377] = precompBuf(K,tau,a,omega,sigma,r,b,c,gammaA,gammaP);
%%%%%
            
%%
ffttime=0;
nltime=0;
tstart=tic
% Dot=zeros(Nr,Nr);
if t==0, 
    erowX=zeros(nl,1); 
    erowY=zeros(nl,1); 
    energyX=zeros(nl,1);
    energyY=zeros(nl,1); 
end;

while round(t/tau)<3e4/tau%6e4
    
    %%% matrix form
%     init(1,:,:)=Epf;
%     init(2,:,:)=Emf;
%     init(3,:,:)=Ppf;
%     init(4,:,:)=Pmf;
%     init(5,:,:)=Npf;
%     init(6,:,:)=Nmf;
%     init(7,:,:)=Mf;
%     ntime=tic;
    [fEpf,fEmf,fPpf,fPmf,fNpf,fNmf,fMf, ep, em]=nonlin2dBuf(Epf,Emf,Ppf,Pmf,Npf,Nmf,Mf,sigma,r);
%     nltime=nltime+toc(ntime);
%     ffttime=ffttime+f1;
    aEpf=exp2M11.*Epf+exp2M12.*Emf+a211.*fEpf+a212.*fEmf;
    aEmf=exp2M22.*Emf+exp2M21.*Epf+a222.*fEmf+a221.*fEpf;
    aPpf=exp2M33.*Ppf+a233.*fPpf;
    aPmf=exp2M44.*Pmf+a244.*fPmf;
	aNpf=exp2M55.*Npf+a255.*fNpf;
	aNmf=exp2M66.*Nmf+a266.*fNmf;
	aMf=exp2M77.*Mf+a277.*fMf;
    
%     ntime=tic;
    [faEpf,faEmf,faPpf,faPmf,faNpf,faNmf,faMf]=nonlin2dBuf(aEpf,aEmf,aPpf,aPmf,aNpf,aNmf,aMf,sigma,r);
%     nltime=nltime+toc(ntime);
%     ffttime=ffttime+f1;
	bEpf=expM11.*Epf+expM12.*Emf+b211.*(2*faEpf-fEpf)+b212.*(2*faEmf-fEmf);
    bEmf=expM22.*Emf+expM21.*Epf+b222.*(2*faEmf-fEmf)+b221.*(2*faEpf-fEpf);
    bPpf=expM33.*Ppf+b233.*(2*faPpf-fPpf);
    bPmf=expM44.*Pmf+b244.*(2*faPmf-fPmf);
	bNpf=expM55.*Npf+b255.*(2*faNpf-fNpf);
	bNmf=expM66.*Nmf+b266.*(2*faNmf-fNmf);
	bMf=expM77.*Mf+b277.*(2*faMf-fMf);
%     ntime=tic;
    [fbEpf,fbEmf,fbPpf,fbPmf,fbNpf,fbNmf,fbMf]=nonlin2dBuf(bEpf,bEmf,bPpf,bPmf,bNpf,bNmf,bMf,sigma,r);
%     nltime=nltime+toc(ntime);
%     ffttime=ffttime+f1;
    cEpf=expM11.*Epf+expM12.*Emf+c111.*fEpf+c112.*fEmf+4*c211.*faEpf+4*c212.*faEmf+c311.*fbEpf+c312.*fbEmf;
    cEmf=expM22.*Emf+expM21.*Epf+c122.*fEmf+c121.*fEpf+4*c222.*faEmf+4*c221.*faEpf+c322.*fbEmf+c321.*fbEpf;
    cPpf=expM33.*Ppf+c133.*fPpf+4*c233.*faPpf+c333.*fbPpf;
    cPmf=expM44.*Pmf+c144.*fPmf+4*c244.*faPmf+c344.*fbPmf;
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
        erowX(nrec)=abs(ep(Nr/2,Nr/2)+em(Nr/2,Nr/2))^2/2;
        erowY(nrec)=abs(1i*(ep(Nr/2,Nr/2)-em(Nr/2,Nr/2)))^2/2;%Ecur(round(Nr/2),round(Nr/2))];
        energyX(nrec)=sum(sum(abs(ep+em).^2));
        energyY(nrec)=sum(sum(abs(1i*(ep-em)).^2));
        if nrec==nl, 
    %         Erow=[Erow erow];erow=[];%DI=[DI dI];dI=[]; % Imax=[Imax imax];imax=[];
    %         disp('write center');
            disp(['t= ',num2str(t),' at ',num2str(toc(tstart)),' s']);
            disp(['sumI= ',num2str((energyX(nrec)+energyY(nrec))/Nr^2)]);
            fid = fopen(datafilenameP, 'a');
            fwrite(fid, erowX, 'double');
            fclose(fid); 
            fid = fopen(datafilenameM, 'a');
            fwrite(fid, erowY, 'double');
            fclose(fid); 
            
            fid = fopen(energyfilenameP, 'a');
            fwrite(fid, energyX, 'double');
            fclose(fid); 
            fid = fopen(energyfilenameM, 'a');
            fwrite(fid, energyY, 'double');
            fclose(fid); 
            nrec=int32(0);   % запись матрицы в файл (40 байт)
            
            if abs(mod(nstep,int32(erow_incr)*10))<tau/2,
%                 figure;colormap(grayCustom);
%                 imagesc(x0,x0,abs(ep+em));colorbar('vert');colormap(grayCustom);%title(['Ix t=',num2str(t)]);
%                 set(gca,'FontSize',gcaFontSize);xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
%                 saveas(gcf,['prof/ex ',filename,' ',num2str(t),'.png'],'png');
% 
%                 imagesc(x0,x0,abs(1i*(ep-em)));colorbar('vert');colormap(grayCustom);%title(['Iy t=',num2str(t)]);
%                 set(gca,'FontSize',gcaFontSize);xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
%                 saveas(gcf,['prof/ey ',filename,' ',num2str(t),'.png'],'png');
% 
%                 imagesc(x0,x0,abs(ep+em) + abs(1i*(ep-em)));colorbar('vert');colormap(grayCustom);%title(['Ix+Iy t=',num2str(t)]);
%                 set(gca,'FontSize',gcaFontSize);xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
%                 saveas(gcf,['prof/e ',filename,' ',num2str(t),'.png'],'png');
% 
%                 imagesc(fftshift(kx),fftshift(kx),abs(fftshift(Epf+Emf)));colorbar('vert');colormap(grayCustom);%title(['F(Ex)t=',num2str(t)]);
%                 set(gca,'FontSize',gcaFontSize);xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
%                 xlim([-10 10]);ylim([-10 10]);saveas(gcf,['prof/exf ',filename,' ',num2str(t),'.png'],'png');
% 
%                 imagesc(fftshift(kx),fftshift(kx),abs(fftshift(1i*(Epf-Emf))));colorbar('vert');colormap(grayCustom);%title(['F(Ey)t=',num2str(t)]);
%                 set(gca,'FontSize',gcaFontSize);xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
%                 xlim([-10 10]);ylim([-10 10]);saveas(gcf,['prof/eyf ',filename,' ',num2str(t),'.png'],'png');
% 
%                 imagesc(fftshift(kx),fftshift(kx),fftshift(abs(Epf+Emf) + abs(1i*(Epf-Emf))));colorbar('vert');colormap(grayCustom);%title(['F(Ex)+F(Ey) t=',num2str(t)]);
%                 set(gca,'FontSize',gcaFontSize);xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
%                 xlim([-10 10]);ylim([-10 10]);saveas(gcf,['prof/ef ',filename,' ',num2str(t),'.png'],'png');
% 
%                 close gcf;
                save(['2D ',filename,'.mat']);
            end
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
    fwrite(fid, erowX(1:nrec), 'double');
    fclose(fid); 
    fid = fopen(datafilenameM, 'a');
    fwrite(fid, erowY(1:nrec), 'double');
    fclose(fid); 
    
    fid = fopen(energyfilenameP, 'a');
    fwrite(fid, energyX(1:nrec), 'double');
    fclose(fid); 
    fid = fopen(energyfilenameM, 'a');
    fwrite(fid, energyY(1:nrec), 'double');
    fclose(fid); 
%     disp(nrec);
    nrec=int32(0);
    erowX=zeros(nl,1);
    erowY=zeros(nl,1);
    energyX=zeros(nl,1);
    energyY=zeros(nl,1);
end;


%%
fid = fopen(datafilenameP, 'rb');
INTp = fread(fid, 'double'); 
fclose(fid); 
subplot('221');
plot((1:numel(INTp))*erow_incr*tau,INTp);
title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
fid = fopen(datafilenameM, 'rb');
INTm = fread(fid, 'double'); 
fclose(fid); 
subplot('222');
plot((1:numel(INTm))*erow_incr*tau,INTm);
fid = fopen(energyfilenameP, 'rb');
energyX = fread(fid, 'double'); 
fclose(fid); 
subplot('223');
plot((1:numel(energyX))*erow_incr*tau,energyX./Nr^2);
fid = fopen(energyfilenameM, 'rb');
energyY = fread(fid, 'double'); 
fclose(fid); 
subplot('224');
plot((1:numel(energyY))*erow_incr*tau,energyY./Nr^2);
gcffilename=['2Dtotal omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' ga=',num2str(gammaA),' gp',num2str(gammaP),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID),'.jpg'];
xlabel('Время t');ylabel('Уровень интенсивности I=|E|^2');saveas(gcf,gcffilename,'jpg');%close gcf;
%%
    figure;imagesc(x0,x0,abs(ep+em));colorbar('vert');colormap(grayCustom);%title(['Ix t=',num2str(t)]);
    set(gca,'FontSize',gcaFontSize);xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
    saveas(gcf,['prof/ex ',filename,' ',num2str(t),'.png'],'png');
    
    figure;imagesc(x0,x0,abs(1i*(ep-em)));colorbar('vert');colormap(grayCustom);%title(['Iy t=',num2str(t)]);
    set(gca,'FontSize',gcaFontSize);xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
    saveas(gcf,['prof/ey ',filename,' ',num2str(t),'.png'],'png');
    
    figure;imagesc(x0,x0,abs(ep+em) + abs(1i*(ep-em)));colorbar('vert');colormap(grayCustom);%title(['Ix+Iy t=',num2str(t)]);
    set(gca,'FontSize',gcaFontSize);xlabel('x','FontSize',18,'FontWeight','bold');ylabel('y','FontSize',18,'FontWeight','bold');
    saveas(gcf,['prof/e ',filename,' ',num2str(t),'.png'],'png');
    
    lim=20;
    figure;imagesc(fftshift(kx),fftshift(kx),abs(fftshift(Epf+Emf)));colorbar('vert');colormap(grayCustom);%title(['F(Ex)t=',num2str(t)]);
    set(gca,'FontSize',gcaFontSize);xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
    xlim([-lim lim]);ylim([-lim lim]);saveas(gcf,['prof/exf ',filename,' ',num2str(t),'.png'],'png');
    
    figure;imagesc(fftshift(kx),fftshift(kx),abs(fftshift(1i*(Epf-Emf))));colorbar('vert');colormap(grayCustom);%title(['F(Ey)t=',num2str(t)]);
    set(gca,'FontSize',gcaFontSize);xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
    xlim([-lim lim]);ylim([-lim lim]);saveas(gcf,['prof/eyf ',filename,' ',num2str(t),'.png'],'png');
    
    figure;imagesc(fftshift(kx),fftshift(kx),fftshift(abs(Epf+Emf) + abs(1i*(Epf-Emf))));colorbar('vert');colormap(grayCustom);%title(['F(Ex)+F(Ey) t=',num2str(t)]);
    set(gca,'FontSize',gcaFontSize);xlabel('k_x','FontSize',18,'FontWeight','bold');ylabel('k_y','FontSize',18,'FontWeight','bold');
    xlim([-lim lim]);ylim([-lim lim]);saveas(gcf,['prof/ef ',filename,' ',num2str(t),'.png'],'png');
%%
dispMode(omega,a,sigma_min,b,c,gammaA,gammaP,r_max);
xf=abs(Epf+Emf);
[m1 x1]=max(xf);
[m2 x2]=max(max(xf));
disp(['kx_max=',num2str(sqrt(kx(x1(x2))^2+kx(x2)^2))]);

yf=abs(1i*(Epf-Emf));
[m1 y1]=max(yf);
[m2 y2]=max(max(yf));
disp(['ky_max=',num2str(sqrt(kx(y1(y2))^2+kx(y2)^2))]);


absX=abs(Epf+Emf).^2;
sumX=sum(sum(absX));
[Xx, KK]=getDistribByWaveNumber2D(absX,kx,256);
figure;plot(KK,Xx/sumX*100,'blue','LineWidth',1,'DisplayName','X')

absY=abs(1i*(Epf-Emf)).^2;
sumY=sum(sum(absY));
[Yy, KK]=getDistribByWaveNumber2D(absY,kx,256);
hold on;grid on;plot(KK,Yy/sumY*100,'red','LineWidth',1,'DisplayName','Y');hold off;
legend('show')
xlim([0 20])%max(KK)
set(gca,'FontSize',gcaFontSize);xlabel('k','FontSize',18,'FontWeight','bold');ylabel('%','FontSize',18,'FontWeight','bold');
saveas(gcf,['prof/NrgByK ',filename,' ',num2str(t),'.png'],'png');
%%
save(['2D ',filename,'.mat']);