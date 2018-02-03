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

% rr=1:0.2:10;
rr=2;
if numel(rr)==1
    dr=1;
else
    dr=rr(2)-rr(1);
end
omega=1;
a=.01;
sigma=0.1;
b=0.01;
c=0.05;
gammaA=-0.05;
gammaP=0.05;

if gammaA ==0 && gammaP==0
    postfix='i';
elseif gammaA ~=0 || gammaP ~= 0
    postfix='a';
else
    postfix='';
end

Ne = 7; % number of equations
omega_aperture_ON = false;
sigma_aperture_ON = false;
r_aperture_ON = false;

k0=real(sqrt(omega/a));
k10=real(sqrt(1/a));
diag=1;
Npolos=100;
if omega > 0
    L=2*pi/k10*Npolos;
else 
    L = 1;
end

h=L/1024;
Nr=ceil(L/h);
stopfilename='stop1';
kfil=int32((Nr-1)/7*3)+1;

r_max=rr(1);sigma_min=sigma;

l=1.2;
d=L;
dx=h;
x0=(0:Nr-1)*h;
X=x0';
% [X,Y]=meshgrid(x0);

% if sigma_aperture_ON==1
%     sigma=sigma_min*1e2*ones(Nr,1);
%     sigma(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=sigma_min;
%     sigma=smooth(sigma,1);
% end

% kx=linspace(-Nr/2,Nr/2,Nr)*2*pi/L;
kx0=((-Nr/2):(Nr/2-1))*2*pi/L;
kx=fftshift(kx0)';
K=kx.^2;
% for i=1:Nr
%     for j=1:Nr
%         K(i,j)=kx(i)^2+kx(j)^2;
%     end
% end
nstep=int32(0);
nrec=int32(0);
erow_incr=10;
nl=1000;
grayCustom=[linspace(1,0,128);linspace(1,0,128);linspace(1,0,128)]';

%%
n0=0.001;
Ep0=(n0*randn(Nr,1));
Em0=(n0*randn(Nr,1));
Pp0=(n0*randn(Nr,1));
Pm0=(n0*randn(Nr,1));
Np0=(n0*randn(Nr,1));
Nm0=(n0*randn(Nr,1));
M0 =(n0*randn(Nr,1));

%%%%%% Single Travelling co-waves X-polarized
% E0=4*b*c*r_max*sigma/(gammaA+sigma)/(b+3*c)-4*b*c*(a*k0^2+gammaP-omega)^2/(b+3*c)/(gammaA+sigma+1)^2-4*b*c/(b+3*c);
% k0=real(sqrt((omega-gammaP)/a));w0=omega+(a*k0^2+gammaP-omega)/(gammaA+sigma+1);
% Ep0=sqrt(E0/2)*exp(1i*(k0*X));
% Em0=sqrt(E0/2)*exp(1i*(k0*X));
% Pp0=sqrt(E0*2)*exp(1i*(k0*X))*(sigma+gammaA-1i*a*k0^2-1i*gammaP-1i*w0)/2/sigma;
% Pm0=sqrt(E0*2)*exp(1i*(k0*X))*(sigma+gammaA-1i*a*k0^2-1i*gammaP-1i*w0)/2/sigma;
% Np0=abs(Ep0).^2*3*(gammaA+sigma)/(4*b*sigma);
% Nm0=abs(Ep0).^2*3*(gammaA+sigma)/(4*b*sigma);
% M0=abs(Ep0).^2*(gammaA+sigma)/(4*c*sigma);

%%%%%% Single Travelling co-waves Y-polarized
% E0=-4*b*c*r_max*sigma/(gammaA-sigma)/(b+3*c)-4*b*c*(-a*k0^2+gammaP+omega)^2/(b+3*c)/(-gammaA+sigma+1)^2-4*b*c/(b+3*c);
% k0=real(sqrt((omega+gammaP)/a));w0=omega-(-a*k0^2+gammaP+omega)/(-gammaA+sigma+1);
% Ep0=-1i*sqrt(E0/2)*exp(1i*(k0*X));
% Em0=1i*sqrt(E0/2)*exp(1i*(k0*X));
% Pp0=-sqrt(E0*2)*exp(1i*(k0*X))*(1i*sigma-1i*gammaA-1i*a*k0^2+gammaP+w0)/2/sigma;
% Pm0=sqrt(E0*2)*exp(1i*(k0*X))*(1i*sigma-1i*gammaA-1i*a*k0^2+gammaP+w0)/2/sigma;
% Np0=abs(Ep0).^2*3*(gammaA-sigma)/(4*b*sigma);
% Nm0=abs(Ep0).^2*3*(gammaA-sigma)/(4*b*sigma);
% M0=abs(Ep0).^2*(gammaA-sigma)/(4*c*sigma);

%%%%%%  Two Travelling trans-waves
% Ep0=sqrt(2*b*c/(3*c+b)*(r-1))*(exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y))));
% Em0=sqrt(2*b*c/(3*c+b)*(r-1))*(exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y))));
% Pp0=Ep0;
% Pm0=Em0;
% Np0=abs(Ep0).^2*3/b*(3+exp(1i*(k0/sqrt(2)*(X+Y)))+exp(-1i*(k0/sqrt(2)*(X+Y)))+0.5*exp(1i*(k0/sqrt(2)*(X-Y)))+0.5*exp(-1i*(k0/sqrt(2)*(X-Y))));
% Nm0=abs(Ep0).^2*3/b*(3+0.5*exp(1i*(k0/sqrt(2)*(X+Y)))+0.5*exp(-1i*(k0/sqrt(2)*(X+Y)))+exp(1i*(k0/sqrt(2)*(X-Y)))+exp(-1i*(k0/sqrt(2)*(X-Y))));
% M0=abs(Ep0).^2/2/c*(exp(1i*(k0/sqrt(2)*(X)))+exp(-1i*(k0/sqrt(2)*(Y)))+exp(1i*(k0/sqrt(2)*(-X)))+exp(-1i*(k0/sqrt(2)*(-Y))));

Epf=fft(Ep0);
Emf=fft(Em0);
Ppf=fft(Pp0);
Pmf=fft(Pm0);
Npf=fft(Np0);
Nmf=fft(Nm0);
Mf =fft(M0);
% Epf(1,1)=Epf(1,1)+1e-4;
Npf=makeReal1d(Npf);
Nmf=makeReal1d(Nmf);
Epf=filter1d(Epf);
Emf=filter1d(Emf);
Ppf=filter1d(Ppf);
Pmf=filter1d(Pmf);
Npf=filter1d(Npf);
Nmf=filter1d(Nmf);
Mf=filter1d(Mf);

mkdir('1d');
mkdir('1d/imag');
mkdir('1d/wc');
shortDirName=['1d/imag/',num2str(gammaA),'/'];
mkdir(shortDirName);
depFileX = ['1d/r-dependence 1d Ex 1..4 omega=',num2str(omega),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP)];
depFileY = ['1d/r-dependence 1d Ey 1..4 omega=',num2str(omega),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP)];
  
%%
while r_max<= max(rr)
    
ID=int16(0);
if round(t/tau) < tau/2
    dirName=[shortDirName,'r',num2str(r_max),postfix];
    mkdir(dirName);
    filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP),' ID#',num2str(ID)];
    while exist(['1d/wc/local 1d x ',filename,'.dat'], 'file')==2,
        ID=ID+1;
        filename=['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' dpi',num2str(Nr),' tau=',num2str(tau),' L=',num2str(L),' ga=',num2str(gammaA),' gp',num2str(gammaP),' ID#',num2str(ID)];
    end
    datafilename=[filename,'.dat'];
    datafilenameP=['1d/wc/local 1d x ',filename,'.dat'];
    fid = fopen(datafilenameP, 'w');
    fclose(fid);
    datafilenameM=['1d/wc/local 1d y ',filename,'.dat'];
    fid = fopen(datafilenameM, 'w');
    fclose(fid);
    energyfilenameP = ['1d/wc/full 1d x ',filename,'.dat'];
    energyfilenameM = ['1d/wc/full 1d y ',filename,'.dat'];
    phasefilenameX = ['1d/wc/phase 1d x ',filename,'.dat'];
    phasefilenameY = ['1d/wc/phase 1d y ',filename,'.dat'];
    refilenameX = ['1d/wc/re 1d x ',filename,'.dat'];
    refilenameY = ['1d/wc/re 1d y ',filename,'.dat'];
end
%%% custom palette

%%% aperture of pumping and losses
% if r_aperture_ON==1
%     r=-1*ones(Nr,Nr);
%     r(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)),round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=r_max;
%     r=smooth(r,1);  
% end
r=r_max;
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
    phaseX=zeros(nl,1);
    phaseY=zeros(nl,1); 
    rerowX=zeros(nl,1); 
    rerowY=zeros(nl,1);
end;
if r_max==rr(1)
    coef=2;
else coef=1;
end
while round(t/tau)*tau<1e4/tau*coef%6e4
    
    %%% matrix form
%     init(1,:,:)=Epf;
%     init(2,:,:)=Emf;
%     init(3,:,:)=Ppf;
%     init(4,:,:)=Pmf;
%     init(5,:,:)=Npf;
%     init(6,:,:)=Nmf;
%     init(7,:,:)=Mf;
%     ntime=tic;
    [fEpf,fEmf,fPpf,fPmf,fNpf,fNmf,fMf, ep, em]=nonlin1d(Epf,Emf,Ppf,Pmf,Npf,Nmf,Mf,sigma,r);
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
    [faEpf,faEmf,faPpf,faPmf,faNpf,faNmf,faMf]=nonlin1d(aEpf,aEmf,aPpf,aPmf,aNpf,aNmf,aMf,sigma,r);
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
    [fbEpf,fbEmf,fbPpf,fbPmf,fbNpf,fbNmf,fbMf]=nonlin1d(bEpf,bEmf,bPpf,bPmf,bNpf,bNmf,bMf,sigma,r);
%     nltime=nltime+toc(ntime);
%     ffttime=ffttime+f1;
    cEpf=expM11.*Epf+expM12.*Emf+c111.*fEpf+c112.*fEmf+4*c211.*faEpf+4*c212.*faEmf+c311.*fbEpf+c312.*fbEmf;
    cEmf=expM22.*Emf+expM21.*Epf+c122.*fEmf+c121.*fEpf+4*c222.*faEmf+4*c221.*faEpf+c322.*fbEmf+c321.*fbEpf;
    cPpf=expM33.*Ppf+c133.*fPpf+4*c233.*faPpf+c333.*fbPpf;
    cPmf=expM44.*Pmf+c144.*fPmf+4*c244.*faPmf+c344.*fbPmf;
	cNpf=expM55.*Npf+c155.*fNpf+4*c255.*faNpf+c355.*fbNpf;
	cNmf=expM66.*Nmf+c166.*fNmf+4*c266.*faNmf+c366.*fbNmf;
	cMf=expM77.*Mf+c177.*fMf+4*c277.*faMf+c377.*fbMf;
    
    Epf=filter1d(cEpf);
    Emf=filter1d(cEmf);
    Ppf=filter1d(cPpf);
    Pmf=filter1d(cPmf);
    Npf=filter1d(cNpf);
    Nmf=filter1d(cNmf);
    Mf=filter1d(cMf);
    
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
        erowX(nrec)=abs(ep(Nr/2)+em(Nr/2))^2/2;
        phaseX(nrec)=angle(ep(Nr/2)+em(Nr/2)/sqrt(2));
        erowY(nrec)=abs(1i*(ep(Nr/2)-em(Nr/2)))^2/2;%Ecur(round(Nr/2),round(Nr/2))];
        phaseY(nrec)=angle(1i*(em(Nr/2)-ep(Nr/2))/sqrt(2));%Ecur(round(Nr/2),round(Nr/2))];
        energyX(nrec)=sum(sum(abs(ep+em).^2))/2;
        energyY(nrec)=sum(sum(abs(1i*(ep-em)).^2))/2;
        rerowX(nrec)=real(ep(Nr/2)+em(Nr/2))/sqrt(2);
        rerowY(nrec)=real(1i*(em(Nr/2)-ep(Nr/2)))/sqrt(2);
        if nrec==nl, 
    %         Erow=[Erow erow];erow=[];%DI=[DI dI];dI=[]; % Imax=[Imax imax];imax=[];
    %         disp('write center');
            disp(['t= ',num2str(t),' at ',num2str(toc(tstart)),' s']);
            disp(['sumI= ',num2str((energyX(nrec)+energyY(nrec))/Nr)]);
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
            
            fid = fopen(phasefilenameX, 'a');
            fwrite(fid, phaseX, 'double');
            fclose(fid); 
            fid = fopen(phasefilenameY, 'a');
            fwrite(fid, phaseY, 'double');
            fclose(fid); 
            
            fid = fopen(refilenameX, 'a');
            fwrite(fid, rerowX, 'double');
            fclose(fid); 
            fid = fopen(refilenameY, 'a');
            fwrite(fid, rerowY, 'double');
            fclose(fid); 
            
            nrec=int32(0);   % ?????? ??????? ? ???? (40 ????)
            
            % Spatial dependence
%             figure;colormap(grayCustom);
%             title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
%             plot(x0,abs(ep+em)/sqrt(2));colorbar('vert');title(['t=',num2str(t)]);
%             saveas(gcf,[dirName,'/ex ',filename,' ',num2str(t),'.png'],'png');
%             plot(x0,abs(1i*(em-ep)/sqrt(2)));colorbar('vert');title(['t=',num2str(t)]);
%             saveas(gcf,[dirName,'/ey ',filename,' ',num2str(t),'.png'],'png');
%             plot(x0,(abs(ep+em) + abs(1i*(em-ep)))/sqrt(2));colorbar('vert');title(['t=',num2str(t)]);
%             saveas(gcf,[dirName,'/e ',filename,' ',num2str(t),'.png'],'png');
%             plot(kx0,abs(fftshift(Epf+Emf)/sqrt(2)));colorbar('vert');title(['t=',num2str(t)]);
%             saveas(gcf,[dirName,'/exf ',filename,' ',num2str(t),'.png'],'png');
%             plot(kx0,abs(fftshift(1i*(Emf-Epf)/sqrt(2))));colorbar('vert');title(['t=',num2str(t)]);
%             saveas(gcf,[dirName,'/eyf ',filename,' ',num2str(t),'.png'],'png');
%             plot(kx0,fftshift(abs(Epf+Emf) + abs(1i*(Emf-Epf)))/sqrt(2));colorbar('vert');title(['t=',num2str(t)]);
%             saveas(gcf,[dirName,'/ef ',filename,' ',num2str(t),'.png'],'png');
%             close gcf;
            
            % Time dependence
%             fid = fopen(datafilenameP, 'rb');
%             INTp = fread(fid, 'double'); 
%             fclose(fid); 
%             subplot('221');
%             plot((1:numel(INTp))*erow_incr*tau,INTp);
%             xlabel('Time t');ylabel('Local I=|Ex|^2');
%             title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
%             fid = fopen(datafilenameM, 'rb');
%             INTm = fread(fid, 'double'); 
%             fclose(fid); 
%             subplot('222');
%             plot((1:numel(INTm))*erow_incr*tau,INTm);
%             xlabel('Time t');ylabel('Local I=|Ey|^2');
%             fid = fopen(energyfilenameP, 'rb');
%             nrgX = fread(fid, 'double'); 
%             fclose(fid); 
%             subplot('223');
%             plot((1:numel(nrgX))*erow_incr*tau,nrgX./Nr);
%             xlabel('Time t');ylabel('Full I=|Ex|^2');
%             fid = fopen(energyfilenameM, 'rb');
%             nrgY = fread(fid, 'double'); 
%             fclose(fid); 
%             subplot('224');
%             plot((1:numel(nrgY))*erow_incr*tau,nrgY./Nr);
%             xlabel('Time t');ylabel('Full I=|Ey|^2');
%             saveas(gcf,[dirName,'/time ',filename,'.png'],'png');
%             close gcf;
% 
            save(['1D ',filename,'.mat']);
            
            if exist(stopfilename,'file')==2,
                break;
            end
        end
%         plot(abs(ep).^2+abs(em).^2);colorbar('vert');
%         drawnow;
    end
    
% % %     plot(1:Nr,abs(ep).^2+abs(em).^2);
% % %     title(t);
% % %     drawnow;
end
toc(tstart);

figure;colormap(grayCustom);
            title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
            plot(x0,abs(ep+em).^2/2);title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/ex ',filename,' ',num2str(t),'.png'],'png');
            plot(x0,abs(1i*(ep-em)).^2/2);title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/ey ',filename,' ',num2str(t),'.png'],'png');
            plot(x0,abs(ep+em).^2/2 + abs(1i*(ep-em)).^2/2);title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/e ',filename,' ',num2str(t),'.png'],'png');
            plot(kx0,abs(fftshift(Epf+Emf)));title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/exf ',filename,' ',num2str(t),'.png'],'png');
            plot(kx0,abs(fftshift(1i*(Epf-Emf))));title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/eyf ',filename,' ',num2str(t),'.png'],'png');
            plot(kx0,fftshift(abs(Epf+Emf) + abs(1i*(Epf-Emf))));title(['t=',num2str(t)]);
            saveas(gcf,[dirName,'/ef ',filename,' ',num2str(t),'.png'],'png');
            close gcf;

            fid = fopen(datafilenameP, 'rb');
            INTp = fread(fid, 'double'); 
            fclose(fid); 
            subplot('321');
            plot((1:numel(INTp))*erow_incr*tau,INTp);
            xlabel('Time t');ylabel('Local I=|Ex|^2');
            title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
            fid = fopen(datafilenameM, 'rb');
            INTm = fread(fid, 'double'); 
            fclose(fid); 
            subplot('322');
            plot((1:numel(INTm))*erow_incr*tau,INTm);
            xlabel('Time t');ylabel('Local I=|Ey|^2');
            fid = fopen(energyfilenameP, 'rb');
            nrgX = fread(fid, 'double'); 
            fclose(fid); 
            subplot('323');
            plot((1:numel(nrgX))*erow_incr*tau,nrgX./Nr);
            xlabel('Time t');ylabel('Full I=|Ex|^2');
            fid = fopen(energyfilenameM, 'rb');
            nrgY = fread(fid, 'double'); 
            fclose(fid); 
            subplot('324');
            plot((1:numel(nrgY))*erow_incr*tau,nrgY./Nr);
            xlabel('Time t');ylabel('Full I=|Ey|^2');
            fid = fopen(phasefilenameX, 'rb');
            phX = fread(fid, 'double'); 
            fclose(fid); 
            subplot('325');
            plot((1:numel(phX))*erow_incr*tau,phX./Nr);
            xlabel('Time t');ylabel('arg(Ex)');
            fid = fopen(phasefilenameY, 'rb');
            phY = fread(fid, 'double'); 
            fclose(fid); 
            subplot('326');
            plot((1:numel(phY))*erow_incr*tau,phY./Nr);
            xlabel('Time t');ylabel('arg(Ey)');
            saveas(gcf,[dirName,'/time ',filename,'.png'],'png');
            close gcf;
            
fid = fopen(depFileX, 'a');
            newNrgX=sum(nrgX(numel(nrgX)-1e4+1:numel(nrgX))./Nr)/1e4;
            fwrite(fid, newNrgX, 'double');
            fclose(fid); 
            fid = fopen(depFileX, 'rb');
            fullX = fread(fid, 'double'); 
            fclose(fid); 
%             plot(rr(1):dr:r_max,fullX,'.-b');
            plot(rr(1)+(0:(numel(fullX)-1))*dr,fullX,'.-b');
            xlabel('Pumping r');ylabel('Full I=|Ex|^2');
            hold on;
            
fid = fopen(depFileY, 'a');
            newNrgY=sum(nrgY(numel(nrgY)-1e4+1:numel(nrgY))./Nr)/1e4;
            fwrite(fid, newNrgY, 'double');
            fclose(fid); 
            fid = fopen(depFileY, 'rb');
            fullY = fread(fid, 'double'); 
            fclose(fid); 
%             plot(rr(1):dr:r_max,fullY,'.-r');
            plot(rr(1)+(0:(numel(fullY)-1))*dr,fullY,'.-r');
            xlabel('Pumping r');ylabel('Local I=|Ey|^2');
            
            saveas(gcf,[shortDirName,'/r',num2str(r_max),' dependence ',filename,'.png'],'png');
            close gcf;
            
r_max=r_max+dr;
t=0;
nrec=int32(0); 

save(['1d/wc/1D ',filename,'.mat']);
            
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
    rerowX=zeros(nl,1);
    rerowY=zeros(nl,1);
end;

%%
save(['1Dbig ',filename,'.mat']);
%%
gcaFontSize=16;
%%
 fid = fopen(datafilenameP, 'rb');
            INTp = fread(fid, 'double'); 
            fclose(fid); 
            subplot('321');
            plot((1:numel(INTp))*erow_incr*tau,INTp);
            xlabel('Time t');ylabel('Local I=|Ex|^2');
            title(['omega=',num2str(omega),' r=',num2str(r_max),' sigma=',num2str(sigma_min),' b=',num2str(b),' a=',num2str(a),' c=',num2str(c),' ga=',num2str(gammaA),' gp',num2str(gammaP)]);
            fid = fopen(datafilenameM, 'rb');
            INTm = fread(fid, 'double'); 
            fclose(fid); 
            subplot('322');
            plot((1:numel(INTm))*erow_incr*tau,INTm);
            xlabel('Time t');ylabel('Local I=|Ey|^2');
            fid = fopen(energyfilenameP, 'rb');
            nrgX = fread(fid, 'double'); 
            fclose(fid); 
            subplot('323');
            plot((1:numel(nrgX))*erow_incr*tau,nrgX./Nr);
            xlabel('Time t');ylabel('Full I=|Ex|^2');
            fid = fopen(energyfilenameM, 'rb');
            nrgY = fread(fid, 'double'); 
            fclose(fid); 
            subplot('324');
            plot((1:numel(nrgY))*erow_incr*tau,nrgY./Nr);
            xlabel('Time t');ylabel('Full I=|Ey|^2');
            fid = fopen(phasefilenameX, 'rb');
            phX = fread(fid, 'double'); 
            fclose(fid); 
            subplot('325');
            plot((1:numel(phX))*erow_incr*tau,phX./Nr);
            xlabel('Time t');ylabel('arg(Ex)');
            fid = fopen(phasefilenameY, 'rb');
            phY = fread(fid, 'double'); 
            fclose(fid); 
            subplot('326');
            plot((1:numel(phY))*erow_incr*tau,phY./Nr);
            xlabel('Time t');ylabel('arg(Ey)');
%% NEAR FIELD X
plot(x0,real((ep+em)/sqrt(2)),'black','LineWidth',1);
xlim([min(x0) max(x0)])
% ylim([0 0.1])
grid on;
set(gca,'FontSize',gcaFontSize);
ylabel('Re(E_x)','FontSize',18,'FontWeight','bold');xlabel('x','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/reEx ',filename,' ',num2str(t),'.png'],'png');

plot(x0,(abs((ep+em)/sqrt(2))).^2,'black','LineWidth',1);
xlim([min(x0) max(x0)])
ylim([0 0.05])
grid on;
set(gca,'FontSize',gcaFontSize);
ylabel('I_x=|E_x|^2','FontSize',18,'FontWeight','bold');xlabel('x','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/ex ',filename,' ',num2str(t),'.png'],'png');
%% FAR FIELD X
plot(fftshift(kx),fftshift(abs((Epf+Emf)/sqrt(2))),'black','LineWidth',1);
xlim([-12 -8])%([min(kx) max(kx)])
% ylim([-0.3 0.3])
grid on;
set(gca,'FontSize',gcaFontSize);
ylabel('F(E_x)','FontSize',18,'FontWeight','bold');xlabel('k_x','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/exf ',filename,' ',num2str(t),'.png'],'png');
%% NEAR FIELD Y
plot(x0,real(1i*(em-ep)/sqrt(2)),'black','LineWidth',1);
xlim([min(x0) max(x0)])
% ylim([0 0.1])
grid on;
set(gca,'FontSize',gcaFontSize);
ylabel('Re(E_y)','FontSize',18,'FontWeight','bold');xlabel('x','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/reEy ',filename,' ',num2str(t),'.png'],'png');

plot(x0,(abs(1i*(em-ep)/sqrt(2))).^2,'black','LineWidth',1);
xlim([min(x0) max(x0)])
ylim([0 0.05])
grid on;
set(gca,'FontSize',gcaFontSize);
ylabel('I_y=|E_y|^2','FontSize',18,'FontWeight','bold');xlabel('x','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/ey ',filename,' ',num2str(t),'.png'],'png');
%% FAR FIELD Y
plot(fftshift(kx),fftshift(abs(1i*(Emf-Epf)/sqrt(2))),'black','LineWidth',1);
xlim([8 12])%xlim([min(kx) max(kx)])
% ylim([-0.3 0.3])
grid on;
set(gca,'FontSize',gcaFontSize);
ylabel('F(E_y)','FontSize',18,'FontWeight','bold');xlabel('k_x','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/eyf ',filename,' ',num2str(t),'.png'],'png');
%%
fid = fopen(refilenameX, 'rb');
reEx=fread(fid, 'double');
fclose(fid); 
fid = fopen(refilenameY, 'rb');
reEy=fread(fid, 'double');
fclose(fid); 
%%
% reEx=real((reEp+reEm))/sqrt(2);
% reEy=real(1i*(reEm-reEp))/sqrt(2);

nw=8192/1;W=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);
plot((numel(reEx)-nw+1:numel(reEx))*erow_incr*tau,reEx(numel(reEx)-nw+1:numel(reEx)),'black');
figure;plot(W,((log(fftshift(abs(fft(reEx(numel(reEx)-nw+1:numel(reEx)))))))),'black')
grid on
xlim([0 min([2 max(W)])]);
set(gca,'FontSize',gcaFontSize);
xlabel('\omega','FontSize',18,'FontWeight','bold');
ylabel('lg(ReEx)','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/spectra Re(Ex) ',filename,'.jpg'],'jpg');%close gcf;
%%
plot((numel(reEy)-nw+1:numel(reEy))*erow_incr*tau,reEy(numel(reEy)-nw+1:numel(reEy)),'black');
figure;plot(W,((log(fftshift(abs(fft(reEy(numel(reEy)-nw+1:numel(reEy)))))))),'black')
xlim([0 min([2 max(W)])]);
set(gca,'FontSize',gcaFontSize);grid on;
xlabel('\omega','FontSize',18,'FontWeight','bold');
ylabel('lg(ReEy)','FontSize',18,'FontWeight','bold');
saveas(gcf,[dirName,'/spectra Re(Ey) ',filename,'.jpg'],'jpg');%close gcf;