%%%% Simple script for VMB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are tuned as infinite and finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON sigma_aperture_ON a K
T=4e5;
t=0;
tau=2e-2;
Nt=ceil(T/tau);

rr=1.3:0.05:4;
dr=rr(2)-rr(1);
omega=1;
a=.01;
sigma=0.1;
b=0.01;
c=0.05;
gammaA=0;
gammaP=0;

if gammaA ==0 && gammaP==0
    postfix='i';
elseif gammaA ~=0 || gammaP ~= 0
    postfix='a';
else
    postfix='';
end

Ne = 7; % number of equations
omega_aperture_ON = false;
sigma_aperture_ON = true;
r_aperture_ON = true;

k0=real(sqrt(omega/a));
k10=real(sqrt(1/a));
diag=1;
Npolos=10;
if omega > 0
    L=2*pi/k10*Npolos;
else 
    L = 1;
end

h=L/256;
Nr=ceil(L/h);
stopfilename='stop1';
kfil=int32((Nr-1)/7*3)+1;

r_max=rr(1);sigma_min=sigma;

l=1.2;
d=L;
dx=h;
x0=(0:Nr-1)*h;
[X,Y]=meshgrid(x0);

if sigma_aperture_ON==1
    sigma=sigma_min*1e2*ones(Nr,Nr);
    sigma(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)),round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=sigma_min;
    sigma=smooth(sigma,1);
end

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
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';

depFileX = ['r-dependence Ex 1..4 omega=',num2str(omega),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP)];
depFileY = ['r-dependence Ey 1..4 omega=',num2str(omega),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP)];
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

%%

while r_max<= max(rr)
    
ID=int16(0);
if t == 0
    dirName=['imag/r',num2str(r_max),postfix];
    mkdir(dirName);
    filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP),' ID#',num2str(ID)];
    while exist(['local2Dx ',filename,'.dat'], 'file')==2,
        ID=ID+1;
        filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP),' ID#',num2str(ID)];
    end
    datafilename=[filename,'.dat'];
    datafilenameP=['local2Dx ',filename,'.dat'];
    fid = fopen(datafilenameP, 'w');
    fclose(fid);
    datafilenameM=['local2Dy ',filename,'.dat'];
    fid = fopen(datafilenameM, 'w');
    fclose(fid);
    energyfilenameP = ['full2Dx ',filename,'.dat'];
    energyfilenameM = ['full2Dy ',filename,'.dat'];
end
%%% custom palette

%%% aperture of pumping and losses
if r_aperture_ON==1
    r=-1*ones(Nr,Nr);
    r(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)),round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=r_max;
    r=smooth(r,1);  
end
%%%%% matrix's precomp
[expM11, expM12, expM21, expM22, expM33, expM44, expM55, expM66, expM77, exp2M11, exp2M12, exp2M21, exp2M22, exp2M33, exp2M44, exp2M55, exp2M66, exp2M77, invM11, invM12, invM21, invM22, invM33, invM44, invM55, invM66, invM77, a211, a212, a221, a222, a233, a244, a255, a266, a277, b211, b212, b221, b222, b233, b244, b255, b266, b277, c111, c112, c121, c122, c133, c144, c155, c166, c177, c211, c212, c221, c222, c233, c244, c255, c266, c277, c311, c312, c321, c322, c333, c344, c355, c366, c377] = precompBuf(K,tau,a,omega,sigma,r,b,c,gammaA,gammaP);
%%%%%

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
if r_max==rr(1)
    coef=30;
else coef=10;
end
while round(t/tau)*tau<1e5*tau*coef%6e4
    
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
        erowX(nrec)=abs(ep(Nr/2,Nr/2)+em(Nr/2,Nr/2))^2;
        erowY(nrec)=abs(1i*(ep(Nr/2,Nr/2)-em(Nr/2,Nr/2)))^2;%Ecur(round(Nr/2),round(Nr/2))];
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
            nrec=int32(0);   % ?????? ??????? ? ???? (40 ????)
            
            % Spatial dependence
            figure;colormap(grayCustom);
            title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
            imagesc(abs(ep+em));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/ex ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(1i*(ep-em)));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/ey ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(ep+em) + abs(1i*(ep-em)));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/e ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(fftshift(Epf+Emf)));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/exf ',filename,' ',num2str(t),'.png'],'png');
            imagesc(abs(fftshift(1i*(Epf-Emf))));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/eyf ',filename,' ',num2str(t),'.png'],'png');
            imagesc(fftshift(abs(Epf+Emf) + abs(1i*(Epf-Emf))));colorbar('vert');title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/ef ',filename,' ',num2str(t),'.png'],'png');
            close gcf;
            
            % Time dependence
            fid = fopen(datafilenameP, 'rb');
            INTp = fread(fid, 'double'); 
            fclose(fid); 
            subplot('221');
            plot((1:numel(INTp))*erow_incr*tau,INTp);
            xlabel('Time t');ylabel('Local I=|Ex|^2');
            title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
            fid = fopen(datafilenameM, 'rb');
            INTm = fread(fid, 'double'); 
            fclose(fid); 
            subplot('222');
            plot((1:numel(INTm))*erow_incr*tau,INTm);
            xlabel('Time t');ylabel('Local I=|Ey|^2');
            fid = fopen(energyfilenameP, 'rb');
            nrgX = fread(fid, 'double'); 
            fclose(fid); 
            subplot('223');
            plot((1:numel(nrgX))*erow_incr*tau,nrgX./Nr^2);
            xlabel('Time t');ylabel('Full I=|Ex|^2');
            fid = fopen(energyfilenameM, 'rb');
            nrgY = fread(fid, 'double'); 
            fclose(fid); 
            subplot('224');
            plot((1:numel(nrgY))*erow_incr*tau,nrgY./Nr^2);
            xlabel('Time t');ylabel('Full I=|Ey|^2');
            saveas(gcf,[dirName,'/time ',filename,'.png'],'png');
            close gcf;

            save(['2D ',filename,'.mat']);
            
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

fid = fopen(depFileX, 'a');
            newNrgX=sum(nrgX(numel(nrgX)-1e4+1:numel(nrgX))./Nr^2)/1e4;
            fwrite(fid, newNrgX, 'double');
            fclose(fid); 
            fid = fopen(depFileX, 'rb');
            fullX = fread(fid, 'double'); 
            fclose(fid); 
            plot(rr(1):dr:r_max,fullX,'.-b');
            xlabel('Pumping r');ylabel('Full I=|Ex|^2');
            hold on;
            
fid = fopen(depFileY, 'a');
            newNrgY=sum(nrgY(numel(nrgY)-1e4+1:numel(nrgY))./Nr^2)/1e4;
            fwrite(fid, newNrgY, 'double');
            fclose(fid); 
            fid = fopen(depFileY, 'rb');
            fullY = fread(fid, 'double'); 
            fclose(fid); 
            plot(rr(1):dr:r_max,fullY,'.-r');
            xlabel('Pumping r');ylabel('Local I=|Ey|^2');
            
            saveas(gcf,['imag/r',num2str(r_max),' dependence ',filename,'.png'],'png');
            close gcf;
            
r_max=r_max+dr;
t=0;
end
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
save(['2Dbig ',filename,'.mat']);