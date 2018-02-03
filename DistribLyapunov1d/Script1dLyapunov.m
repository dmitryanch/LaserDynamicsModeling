%%%% Simple script for MB-system by etdrk3-method. Sigma, Delta and R 
%%%% parameters are tuned as infinite and finite aperture

clc
clear
%%
global x0 dx x Nr kfil r_aperture_ON delta_aperture_ON sigma_aperture_ON a K
T=4e5;
t=0;
tau=1e-3;
Nt=ceil(T/tau);
sigma=1;
gamma=2.2;
delta=0;
a=.01;
r=168.75;

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
    r=-1*ones(Nr,1);
    r(round(Nr/l*(l/2-0.5)):round(Nr/l*(l/2+0.5)))=r_max;
    r=smooth(r,1);  
end
if sigma_aperture_ON==1
    sigma=sigma_min*5*ones(Nr,1);
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
E0=1/sqrt(2)*ones(Nr,1).*sqrt(r_max-1-(delta./(1+sigma_min)).^2).*r/r_max+1i/sqrt(2)*ones(Nr,1).*sqrt(r_max-1-(delta./(1+sigma_min)).^2).*r/r_max;
P0=E0*(1-1i*delta./(1+sigma_min));
D0=ones(Nr,1).*(1+(delta./(1+sigma_min)).^2).*r/r_max;
E0(1,1)=E0(1,1)+1e-4;

Pf=fft(P0);
Ef=fft(E0);
Df=fft(D0);
Ef(1,1)=Ef(1,1)+1e-4;
Df=makeReal1d(Df);
Ef=filter1d(Ef);
Pf=filter1d(Pf);
Df=filter1d(Df);

Ef3=Ef;
Pf3=Pf;
Df3=Df;

% kx=linspace(-Nr/2,Nr/2,Nr)*2*pi/L;
kx=((-Nr/2):(Nr/2-1))*2*pi/L;
kx=fftshift(kx);
K=zeros(Nr,1);
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
    cc1(:,1)=(-4-ccc*tau+exp(ccc*tau).*(4-3*ccc*tau+(ccc*tau).^2))./(ccc.^3*tau^2);
    cc1(:,2)=(2+ccc*tau+exp(ccc*tau).*(-2+ccc*tau))./(ccc.^3*tau^2);
    cc1(:,3)=(-4-3*ccc*tau-(ccc*tau).^2+exp(ccc*tau).*(4-ccc*tau))./(ccc.^3*tau^2);

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
cc=zeros(5,5,Nr);ecc=zeros(5,5,Nr);ecc2=zeros(5,5,Nr);
for q=1:Nr
    cc(:,:,q)=[-sigma,a*K(q),sigma,0,0;
        -a*K(q),-sigma,0,sigma,0;
        0,0,-1,delta,0;
        0,0,-delta,-1,0;
        0,0,0,0,-gamma];
    if sigma_aperture_ON==1
        ccclr = -1;
        cc(1,1,q)=0;
        cc(1,3,q)=0;
        cc(2,2,q)=0;
        cc(2,4,q)=0;
    end
    ecc(:,:,q)=expm(cc(:,:,q)*tau*2);
    ecc2(:,:,q)=expm(cc(:,:,q)*tau);
end
%     invc=inv(c);
    
    
    ttau=2*tau;
    
    ccc1=zeros(5,5,Nr);ccc2=zeros(5,5,Nr);ccc3=zeros(5,5,Nr);
    for q=1:Nr
%         ccc1(:,:,q)=(-4*eye(5)-ttau*cc(:,:,q)+ecc(:,:,q)*(4*eye(5)-3*ttau*cc(:,:,q)+ttau^2*cc(:,:,q)^2))/ttau^2/(cc(:,:,q)^3);
%         ccc2(:,:,q)=(2*eye(5)+ttau*cc(:,:,q)+ecc(:,:,q)*(-2*eye(5)+ttau*cc(:,:,q)))/ttau^2/(cc(:,:,q)^3);
%         ccc3(:,:,q)=(-4*eye(5)-3*ttau*cc(:,:,q)-ttau^2*cc(:,:,q)^2+ecc(:,:,q)*(4*eye(5)-ttau*cc(:,:,q)))/ttau^2/(cc(:,:,q)^3);
        
        ccc1(:,:,q) = -4/ttau^2*(cc(:,:,q)^-1)^3 - 1/ttau*(cc(:,:,q)^-1)^2 + ecc(:,:,q)*4/ttau^2*(cc(:,:,q)^-1)^3 - ecc(:,:,q)*3/ttau*(cc(:,:,q)^-1)^2 + ecc(:,:,q)*(cc(:,:,q)^-1);
        ccc2(:,:,q) = 2/ttau^2*(cc(:,:,q)^-1)^3 + 1/ttau*(cc(:,:,q)^-1)^2 + ecc(:,:,q)*(-2)/ttau^2*(cc(:,:,q)^-1)^3 + ecc(:,:,q)/ttau*(cc(:,:,q)^-1)^2;
        ccc3(:,:,q) = -4/ttau^2*(cc(:,:,q)^-1)^3 - 3/ttau*(cc(:,:,q)^-1)^2 - cc(:,:,q)^-1 + ecc(:,:,q)*4/ttau^2*(cc(:,:,q)^-1)^3 - ecc(:,:,q)/ttau*(cc(:,:,q)^-1)^2;
    end
    
%     cc1L(1)=(-4*eye(5)-ccclr*ttau+exp(ccclr*ttau).*(4-3*ccclr*ttau+(ccclr*ttau).^2))./(ccclr.^3*ttau^2);
%     cc1L(2)=(2+ccclr*ttau+exp(ccclr*ttau).*(-2+ccclr*ttau))./(ccclr.^3*ttau^2);
%     cc1L(3)=(-4-3*ccclr*ttau-(ccclr*ttau).^2+exp(ccclr*ttau).*(4-ccclr*ttau))./(ccclr.^3*ttau^2);
% 
%     if sigma_aperture_ON==1
%         cccli = -1;
%     else
%         cccli = -sigma;
%     end
%     cc2L(1)=(-4-cccli*ttau+exp(cccli*ttau).*(4-3*cccli*ttau+(cccli*ttau).^2))./(cccli.^3*ttau^2);
%     cc2L(2)=(2+cccli*ttau+exp(cccli*ttau).*(-2+cccli*ttau))./(cccli.^3*ttau^2);
%     cc2L(3)=(-4-3*cccli*ttau-(cccli*ttau).^2+exp(cccli*ttau).*(4-cccli*ttau))./(cccli.^3*ttau^2);
% 
%     dddlr=-1.0;
%     cc3L(1)    =(-4-dddlr*ttau+exp(dddlr*ttau).*(4-3*dddlr*ttau+(dddlr*ttau).^2))./(dddlr.^3*ttau^2);
%     cc3L(2)    =(2+dddlr*ttau+exp(dddlr*ttau).*(-2+dddlr*ttau))./(dddlr.^3*ttau^2);
%     cc3L(3)    =(-4-3*dddlr*ttau-(dddlr*ttau).^2+exp(dddlr*ttau).*(4-dddlr*ttau))./(dddlr.^3*ttau^2);
%   
%     dddli=-1.0;
%     cc4L(1)    =(-4-dddli*ttau+exp(dddli*ttau).*(4-3*dddli*ttau+(dddli*ttau).^2))./(dddli.^3*ttau^2);
%     cc4L(2)    =(2+dddli*ttau+exp(dddli*ttau).*(-2+dddli*ttau))./(dddli.^3*ttau^2);
%     cc4L(3)    =(-4-3*dddli*ttau-(dddli*ttau).^2+exp(dddli*ttau).*(4-dddli*ttau))./(dddli.^3*ttau^2);
%   
%     eeel=-gamma;
%     cc5L(1)	  =(-4-eeel*ttau+exp(eeel*ttau).*(4-3*eeel*ttau+(eeel*ttau).^2))./(eeel.^3*ttau^2);
%     cc5L(2)    =(2+eeel*ttau+exp(eeel*ttau).*(-2+eeel*ttau))./(eeel.^3*ttau.^2);
%     cc5L(3)    =(-4-3*eeel*ttau-(eeel*ttau).^2+exp(eeel*ttau).*(4-eeel*ttau))./(eeel.^3*ttau^2);

    e0r = ifft(filter1d(rand(1,Nr)));
    e0i = ifft(filter1d(rand(1,Nr)));
    p0r = ifft(filter1d(rand(1,Nr)));
    p0i = ifft(filter1d(rand(1,Nr)));
    d0r = ifft(filter1d(rand(1,Nr)));

%     norm_0=norm([e0r; e0i; p0r; p0i; d0r],'inf');
    norm_0=norm([e0r, e0i, p0r, p0i, d0r],'inf');
    e0r=e0r/norm_0(1); e0i=e0i/norm_0(1); p0r=p0r/norm_0(1); p0i=p0i/norm_0(1); d0r=d0r/norm_0(1);
    
    p_r=fft(p0r);
    p_i=fft(p0i);
    e_r=fft(e0r);
    e_i=fft(e0i);
    d_r=fft(d0r);
    
%%%%%
A=zeros(5,Nr);
B=zeros(5,Nr);
C=zeros(5,Nr);
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
    
    Ef=Ef.*exp(ccc*tau)+cc1(:,1).*FEf +4*FAEf.*cc1(:,2)+FBEf.*cc1(:,3);
    Pf=Pf.*exp(ddd*tau)+cc2(1).*FPf     +4*FAPf.*cc2(2)    +FBPf.*cc2(3);
    Df=Df.*exp(eee*tau)+cc3(1).*FDf     +4*FADf.*cc3(2)    +FBDf.*cc3(3);
    
    Ef=filter1d(Ef);
    Pf=filter1d(Pf);
    Df=filter1d(Df);
    
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
            
%             %%% near field images
%             imagesc(abs(ifft(Ef)));colorbar('vert');colormap(gray);
%             gcffilename=[filename,' nf ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
%             colormap(grayCustom);
%             gcffilename=[filename,' nf inv ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
%             %%% far field images
%             imagesc(abs(fftshift(Ef)));colorbar('vert');colormap(gray);
%             gcffilename=[filename,' ff ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
%             colormap(grayCustom);
%             gcffilename=[filename,' ff inv ',num2str(round(t/tau)),'.png'];
%             saveas(gcf,gcffilename,'png');
%             %%% closing gcf
% %             close gcf;
        end
    end
    
    %%%%% Lyapunov's exponents calculating
    
    if lyapunov_ON == 1 && mod(nstep,2)==0 && nstep ~= 0
         
        init=[e_r;e_i;p_r;p_i;d_r];
        [Fe_r, Fe_i, Fp_r, Fp_i, Fd_r]=nonlinear1dLyapunov(e_r, e_i, p_r, p_i, d_r, Ef3,Pf3,Df3,sigma,delta,gamma,r);       
        F=[Fe_r; Fe_i; Fp_r; Fp_i; Fd_r];
        
        for q=1:Nr
            A(:,q)=ecc2(:,:,q)*init(:,q)+(ecc2(:,:,q)-eye(5))/cc(:,:,q)*F(:,q);
        end
        
%         Aefr=e_r.*exp(0.5*ccclr*ttau)+(exp(0.5*ccclr*ttau)-1).*Fe_r./ccclr;
%         Aefi=e_i.*exp(0.5*cccli*ttau)+(exp(0.5*cccli*ttau)-1).*Fe_i./cccli;
%         Apfr=p_r.*exp(0.5*dddlr*ttau)+(exp(0.5*dddlr*ttau)-1).*Fp_r./dddlr;
%         Apfi=p_i.*exp(0.5*dddli*ttau)+(exp(0.5*dddli*ttau)-1).*Fp_i./dddli;
%         Adf=d_r.*exp(0.5*eeel*ttau)+(exp(0.5*eeel*ttau)-1).*Fd_r./eeel;

%         [FAe_r, FAe_i, FAp_r, FAp_i, FAd_r]=nonlinear1dLyapunov(Aefr, Aefi, Apfr, Apfi, Adf, Ef2,Pf2,Df2,sigma,delta,gamma,r);
        [FAe_r, FAe_i, FAp_r, FAp_i, FAd_r]=nonlinear1dLyapunov(A(1,:), A(2,:), A(3,:), A(4,:), A(5,:), Ef2,Pf2,Df2,sigma,delta,gamma,r);
        FA=[FAe_r; FAe_i; FAp_r; FAp_i; FAd_r];
        for q=1:Nr
            B(:,q)=ecc(:,:,q)*init(:,q)+(ecc(:,:,q)-eye(5))/cc(:,:,q)*(2*FA(:,q)-F(:,q));
        end
        
%         Befr=e_r.*exp(ccclr*ttau)+(exp(ccclr*ttau)-1).*(2*FAe_r-Fe_r)./ccclr;
%         Befi=e_i.*exp(cccli*ttau)+(exp(cccli*ttau)-1).*(2*FAe_i-Fe_i)./cccli;
%         Bpfr=p_r.*exp(dddlr*ttau)+(exp(dddlr*ttau)-1).*(2*FAp_r-Fp_r)./dddlr;
%         Bpfi=p_i.*exp(dddli*ttau)+(exp(dddli*ttau)-1).*(2*FAp_i-Fp_i)./dddli;
%         Bdf=d_r.*exp(eeel*ttau)+(exp(eeel*ttau)-1).*(2*FAd_r-Fd_r)./eeel;

        [FBe_r, FBe_i, FBp_r, FBp_i, FBd_r]=nonlinear1dLyapunov(B(1,:), B(2,:), B(3,:), B(4,:), B(5,:), Ef,Pf,Df,sigma,delta,gamma,r);
        FB=[FBe_r; FBe_i; FBp_r; FBp_i; FBd_r];
        for q=1:Nr
            C(:,q)=ecc(:,:,q)*init(:,q)+ccc1(:,:,q)*F(:,q)+4*ccc2(:,:,q)*FA(:,q)+ccc3(:,:,q)*FB(:,q);
        end
        
%         e_r=e_r.*exp(ccclr*ttau)  +cc1L(1).*Fe_r   +4*FAe_r.*cc1L(2)  +FBe_r.*cc1L(3);
%         e_i=e_i.*exp(cccli*ttau)  +cc2L(1).*Fe_i   +4*FAe_i.*cc2L(2)  +FBe_i.*cc2L(3);
%         p_r=p_r.*exp(dddlr*ttau)  +cc3L(1).*Fp_r   +4*FAp_r.*cc3L(2)  +FBp_r.*cc3L(3);
%         p_i=p_i.*exp(dddli*ttau)  +cc4L(1).*Fp_i   +4*FAp_i.*cc4L(2)  +FBp_i.*cc4L(3);
%         d_r=d_r.*exp(eeel*ttau)   +cc5L(1).*Fd_r   +4*FAd_r.*cc5L(2)  +FBd_r.*cc5L(3);
        
        e_r=filter1d(C(1,:));
        e_i=filter1d(C(2,:));
        p_r=filter1d(C(3,:));
        p_i=filter1d(C(4,:));
        d_r=filter1d(C(5,:));
        
        
        if abs(mod(nstep,int32(10)))<tau/2,
            er=ifft(e_r);
            ei=ifft(e_i);
            pr=ifft(p_r);
            pi=ifft(p_i);
            dr=ifft(d_r);
            
%             norm_vari=norm([er; ei; pr; pi; dr],'inf');
            norm_vari=norm([er, ei, pr, pi, dr],'inf');
            PP = log ( norm_vari );
            SP = SP + PP;
            
            e_r = fft(er / norm_vari); 
            e_i = fft(ei / norm_vari); 
            p_r = fft(pr / norm_vari); 
            p_i = fft(pi / norm_vari); 
            d_r = fft(dr / norm_vari);        
    
            if abs(mod(nstep,int32(2000)))<tau/2,
                disp(['t= ',num2str(t),' SP= ',num2str(SP / t)]);
                SPar=[SPar SP / t];
                subplot(211);
                plot(linspace(0,t,numel(SPar)),SPar);
                if nstep > nl*erow_incr
                    subplot(212);
                    INT=[];
                    fid = fopen(datafilename, 'rb');
                    INT = fread(fid, 'double'); 
                    fclose(fid); 
                    % Ist=r_max-1-(delta/(1+sigma_min))^2;
                    plot((1:numel(INT))*erow_incr*tau,INT);
                end
                drawnow;
            end
        end
        Ef3 = Ef;
        Pf3 = Pf;
        Df3 = Df;
    else
        Ef2 = Ef;
        Pf2 = Pf;
        Df2 = Df;
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
