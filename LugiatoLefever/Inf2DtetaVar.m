%%%% Simple script for MB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are can be set as on infinite or on finite aperture

clc
clear
%%
global x0 dx x Nr kfil
T=4e5;
t=0;
tau=0.005;
Nt=ceil(T/tau);
sigma=2;
b2=0.1;  % diffraction constant. Its square
Ein=2; % real and positive >=0
sigma0=sigma;
L=1;
h=L/128; 	% Spatial resolution
Nr=ceil(L/h);
stopfilename='stop0';        %%% file name that will stop calculating if exist
kfil=int32((Nr-1)/7*3)+1;     %%% parameter for dealising

kx=((-Nr/2):(Nr/2-1))*2*pi/L;	
kx=fftshift(kx);
ky=kx;										
% ky=fftshift(ky);
K=zeros(Nr,Nr);
for i=1:Nr
    for j=1:Nr
        K(i,j)=kx(i)^2+ky(j)^2;		%%% array of wave numbers
    end
end

nstep=int32(0);
nrec=int32(0);
erow_incr=10;
nl=1000;

d=L;
dx=h;
% x0=linspace(-d/2,d/2,Nr);
x0=(0:Nr-1)*h;
[X Y]=meshgrid(x0);

%%
E0=0.1*rand(Nr,Nr);

Ef=fft2(E0);			%%%  instant distibution of Furier-coefficient
Ef(1,1)=Ef(1,1)+1e-4;

Ef=filter2d(Ef);
%%
ID=int16(0);
mkdir('time');			%%% some directories creates for saving pictures and numerics data (intensity, real part of field, total intesity by whole surface)
mkdir('time/2d');
picsDir='pics/2d';
mkdir(picsDir);
filename=[' sigma=',num2str(sigma),' Ein=',num2str(Ein),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
while exist(['time/2d/INT ',filename,'.dat'], 'file')==2,
    ID=ID+1;
    filename=[' sigma=',num2str(sigma),' Ein=',num2str(Ein),' dpi',num2str(Nr),' tau=',num2str(tau),' t=',num2str(t),' L=',num2str(L),' ID#',num2str(ID)];
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

%%%%% precompilation extra expressions for performance improvement
    
c=-1 + 1i*b2*K;
hc=c*tau;
exp2M1=exp(hc/2);
exp2M2=(exp(hc/2)-1)./c;
expM1=exp(hc);
expM2=(exp(hc)-1)./c;
c1=(-4./c.^3-tau./c.^2+exp(hc).*(4./c.^3-3*tau./c.^2+tau^2./c))/tau/tau;
c2=(2./c.^3+tau./c.^2+exp(hc).*(-2./c.^3+tau./c.^2))/tau/tau;
c3=(-4./c.^3-3*tau./c.^2-tau^2./c+exp(hc).*(4./c.^3-tau./c.^2))/tau/tau;

%%%%% calculating loop
while round(t/tau)<10e4/tau			%%% calculating time restriction HERE
    sigma = sigma0;     %%% here you should set sigma(t)
    
    [FEf,e1]=nonlinear2dInfTetaVar(Ef,Ein,sigma);			%%% nonlinear part of equation
    AEf = exp2M1.*Ef + exp2M2.*FEf;				
     
    [FAEf]=nonlinear2dInfTetaVar(AEf,Ein,sigma);
    BEf = expM1.*Ef + expM2.*(2*FAEf-FEf);
    
    [FBEf]=nonlinear2dInfTetaVar(BEf,Ein,sigma);
    cEf=expM1.*Ef + c1.*FEf + 4*c2.*FAEf + c3.*FBEf;
    
    Ef=filter2d(cEf);						%%% dealising
    
    nstep=nstep+int32(1);
    t=t+tau;
	
	%%% saving real part of field, intesity and total intensity (YOU CAN COMMENT THIS PART IF YOU DONT NEED IT)
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
	
	
	%%% here you can choose something to display
%     imagesc(abs(e1).^2);colorbar('vert');				%%% distribution of intensity in physical dimensions (NEAR FIELD)
% % %     imagesc(abs(fftshift(Ef)).^2);colorbar('vert');		%%% distribution of intensity of Furier-components (FAR FIELD)
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
saveas(gcf,gcffilename,'jpg');%close gcf;

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
%%  FREQUENCY DEPENDENCE OF REAL PART OF FIELD
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
% saveas(gcf,['NrgByK ',filename,' ',num2str(t),'.png'],'png');

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

%% FREQUENCY DEPENDENCE OF INTESITY IN RANDOM POINT
nw=min(8192/4,numel(INT));W=(-nw/2:nw/2-1)*2*pi/(tau*nw*erow_incr);
plot((numel(INT)-nw+1:numel(INT))*erow_incr*tau,INT(numel(INT)-nw+1:numel(INT)),'black');
figure;plot(W,log((fftshift(abs(fft(INT(numel(INT)-nw+1:numel(INT)))))).^0.1),'black')
xlim([0 min([1.5 max(W)])]);
xlabel('\omega','FontSize',18,'FontWeight','bold');
ylabel('lg(C)','FontSize',18,'FontWeight','bold');
% saveas(gcf,[picsDir,'/spectra  ',filename,'.jpg'],'jpg');%close gcf;