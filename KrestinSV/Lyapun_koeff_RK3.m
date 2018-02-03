clc
clear
%%
global E P D
global e_r e_i p_r p_i d_r
global SP
T=2e6;
tau=1;
tau_2=tau/2;
Nt=ceil(T/tau);
Nt_beg=1;
Nt_end=1e7;
L=2*pi/35.5884443058738*12;
h=L/128;
Nr=ceil(L/h);
gamma=2.2;
sigma=1;
delta=0;
a=0.01;
r=60;

disp(Nr);
disp(Nt);

% sigma=5
% gamma=0.1
% delta=-3
% r=30
% a = 0.01
%%

E0=sqrt(r-1-(delta/(1+sigma))^2)*ones(Nr,1);E0(1,1)=1.1*E0(1,1); %E0(1,1)=E0(1,1)+1e-2;
P0=E0*(1-1i*delta/(1+sigma));
D0=(1+(delta/(1+sigma))^2)*ones(Nr,1);
D0(1,1) = 1.1 * D0(1,1);

% E0=0.001*rand(Nr,1)+0.001*1i*rand(Nr,1);
% P0=0.001*rand(Nr,1)+0.001*1i*rand(Nr,1);
% D0=0.001*rand(Nr,1);

e_r0=rand(Nr,1);
e_i0=rand(Nr,1);
p_r0=rand(Nr,1);
p_i0=rand(Nr,1);
d_r0=rand(Nr,1);

norm_0=norm([e_r0; e_i0; p_r0; p_i0; d_r0],'inf');

e_r0=e_r0/norm_0; e_i0=e_i0/norm_0; p_r0=p_r0/norm_0; p_i0=p_i0/norm_0; d_r0=d_r0/norm_0;

disp(norm([e_r0; e_i0; p_r0; p_i0; d_r0],'inf'));

P=P0;
E=E0;
D=D0;

e_r=e_r0;
e_i=e_i0;
p_r=p_r0;
p_i=p_i0;
d_r=d_r0;

kx1=(((-Nr/2):(Nr/2-1))*2*pi/L)';
kx=fftshift(kx1);
K=kx.^2;
%%
tic;


SP = 0;

% Nt_beg=1e7+1;
% Nt_end=2e7;

for m=Nt_beg:Nt_end
    
    E_f=fft(E); P_f=fft(P); D_f=fft(D);
    
    e_r_f=fft(e_r); e_i_f=fft(e_i); p_r_f=fft(p_r); p_i_f=fft(p_i); d_r_f=fft(d_r);
    
    f_er=sigma*(p_r_f-e_r_f)+a*K.*e_i_f;
    
    f_ei=sigma*(p_i_f-e_i_f)-a*K.*e_r_f;
    
    f_pr=delta*p_i_f-p_r_f+fft(D.*e_r+real(E).*d_r);
    
    f_pi=fft(D.*e_i+imag(E).*d_r)-delta*p_r_f-p_i_f;
    
    f_dr=-gamma*fft(real(P).*e_r+imag(P).*e_i+real(E).*p_r+imag(E).*p_i+d_r);
    
    fi0_er=tau*f_er; fi0_ei=tau*f_ei; fi0_pr=tau*f_pr; fi0_pi=tau*f_pi; fi0_dr=tau*f_dr;

    F_E=-1i*a*K.*E_f+sigma*(P_f-E_f);

    DE=D.*E;

    F_P=-(1+1i*delta)*P_f+fft(DE);

    EP=(conj(E).*P+E.*conj(P))/2;

    F_D=-gamma*fft(D-r+EP);

    Fi0_E=tau_2*F_E; Fi0_P=tau_2*F_P; Fi0_D=tau_2*F_D;

    Yt_E=E_f+Fi0_E/2; Yt_P=P_f+Fi0_P/2; Yt_D=D_f+Fi0_D/2; 

    F_E=-1i*a*K.*Yt_E+sigma*(Yt_P-Yt_E);

    DE=ifft(Yt_D).*ifft(Yt_E);

    F_P=-(1+1i*delta)*Yt_P+fft(DE);

    iff_E=ifft(Yt_E);

    iff_P=ifft(Yt_P);

    EP=(conj(iff_E).*iff_P+iff_E.*conj(iff_P))/2;

    F_D=-gamma*(Yt_D+fft(EP-r));

    Fi1_E=tau_2*F_E; Fi1_P=tau_2*F_P; Fi1_D=tau_2*F_D;

    Yt_E=E_f-Fi0_E+2*Fi1_E; Yt_P=P_f-Fi0_P+2*Fi1_P; Yt_D=D_f-Fi0_D+2*Fi1_D;

    F_E=-1i*a*K.*Yt_E+sigma*(Yt_P-Yt_E);

    DE=ifft(Yt_D).*ifft(Yt_E);

    F_P=-(1+1i*delta)*Yt_P+fft(DE);

    iff_E=ifft(Yt_E);

    iff_P=ifft(Yt_P);

    EP=(conj(iff_E).*iff_P+iff_E.*conj(iff_P))/2;

    F_D=-gamma*(Yt_D+fft(EP-r));

    Fi2_E=tau_2*F_E; Fi2_P=tau_2*F_P; Fi2_D=tau_2*F_D;

    E_f=E_f+(Fi0_E+4*Fi1_E+Fi2_E)/6;

    P_f=P_f+(Fi0_P+4*Fi1_P+Fi2_P)/6;

    D_f=D_f+(Fi0_D+4*Fi1_D+Fi2_D)/6;

    E=ifft(E_f); P=ifft(P_f); D=ifft(D_f);
    
    Yt_er=e_r_f+fi0_er/2; Yt_ei=e_i_f+fi0_ei/2; Yt_pr=p_r_f+fi0_pr/2; Yt_pi=p_i_f+fi0_pi/2; Yt_dr=d_r_f+fi0_dr/2;
        
    f_er=sigma*(Yt_pr-Yt_er)+a*K.*Yt_ei;
    
    f_ei=sigma*(Yt_pi-Yt_ei)-a*K.*Yt_er;    
    
    iff_er=ifft(Yt_er);
    
    iff_ei=ifft(Yt_ei);
    
    iff_dr=ifft(Yt_dr);
    
    f_pr=delta*Yt_pi-Yt_pr+fft(D.*iff_er+real(E).*iff_dr);
    
    f_pi=fft(D.*iff_ei+imag(E).*iff_dr)-delta*Yt_pr-Yt_pi;        
    
    f_dr=-gamma*Yt_dr-gamma*fft(real(P).*iff_er+imag(P).*iff_ei+real(E).*ifft(Yt_pr)+imag(E).*ifft(Yt_pi));
    
    fi1_er=tau*f_er; fi1_ei=tau*f_ei; fi1_pr=tau*f_pr; fi1_pi=tau*f_pi; fi1_dr=tau*f_dr;
    
    F_E=-1i*a*K.*E_f+sigma*(P_f-E_f);

    DE=D.*E;

    F_P=-(1+1i*delta)*P_f+fft(DE);

    EP=(conj(E).*P+E.*conj(P))/2;

    F_D=-gamma*fft(D-r+EP);

    Fi0_E=tau_2*F_E; Fi0_P=tau_2*F_P; Fi0_D=tau_2*F_D;

    Yt_E=E_f+Fi0_E/2; Yt_P=P_f+Fi0_P/2; Yt_D=D_f+Fi0_D/2; 

    F_E=-1i*a*K.*Yt_E+sigma*(Yt_P-Yt_E);

    DE=ifft(Yt_D).*ifft(Yt_E);

    F_P=-(1+1i*delta)*Yt_P+fft(DE);

    iff_E=ifft(Yt_E);

    iff_P=ifft(Yt_P);

    EP=(conj(iff_E).*iff_P+iff_E.*conj(iff_P))/2;

    F_D=-gamma*(Yt_D+fft(EP-r));

    Fi1_E=tau_2*F_E; Fi1_P=tau_2*F_P; Fi1_D=tau_2*F_D;

    Yt_E=E_f-Fi0_E+2*Fi1_E; Yt_P=P_f-Fi0_P+2*Fi1_P; Yt_D=D_f-Fi0_D+2*Fi1_D;

    F_E=-1i*a*K.*Yt_E+sigma*(Yt_P-Yt_E);

    DE=ifft(Yt_D).*ifft(Yt_E);

    F_P=-(1+1i*delta)*Yt_P+fft(DE);

    iff_E=ifft(Yt_E);

    iff_P=ifft(Yt_P);

    EP=(conj(iff_E).*iff_P+iff_E.*conj(iff_P))/2;

    F_D=-gamma*(Yt_D+fft(EP-r));

    Fi2_E=tau_2*F_E; Fi2_P=tau_2*F_P; Fi2_D=tau_2*F_D;

    E_f=E_f+(Fi0_E+4*Fi1_E+Fi2_E)/6;

    P_f=P_f+(Fi0_P+4*Fi1_P+Fi2_P)/6;

    D_f=D_f+(Fi0_D+4*Fi1_D+Fi2_D)/6;

    E=ifft(E_f); P=ifft(P_f); D=ifft(D_f);
    
    Yt_er=e_r_f-fi0_er+2*fi1_er; Yt_ei=e_i_f-fi0_ei+2*fi1_ei; Yt_pr=p_r_f-fi0_pr+2*fi1_pr; Yt_pi=p_i_f-fi0_pi+2*fi1_pi; Yt_dr=d_r_f-fi0_dr+2*fi1_dr;
        
    f_er=sigma*(Yt_pr-Yt_er)+a*K.*Yt_ei;
    
    f_ei=sigma*(Yt_pi-Yt_ei)-a*K.*Yt_er;    
    
    iff_er=ifft(Yt_er);
    
    iff_ei=ifft(Yt_ei);
    
    iff_dr=ifft(Yt_dr);
    
    f_pr=delta*Yt_pi-Yt_pr+fft(D.*iff_er+real(E).*iff_dr);
    
    f_pi=fft(D.*iff_ei+imag(E).*iff_dr)-delta*Yt_pr-Yt_pi;        
    
    f_dr=-gamma*Yt_dr-gamma*fft(real(P).*iff_er+imag(P).*iff_ei+real(E).*ifft(Yt_pr)+imag(E).*ifft(Yt_pi));
    
    fi2_er=tau*f_er; fi2_ei=tau*f_ei; fi2_pr=tau*f_pr; fi2_pi=tau*f_pi; fi2_dr=tau*f_dr;
    
    e_r_f=e_r_f+(fi0_er+4*fi1_er+fi2_er)/6;
    
    e_i_f=e_i_f+(fi0_ei+4*fi1_ei+fi2_ei)/6;
    
    p_r_f=p_r_f+(fi0_pr+4*fi1_pr+fi2_pr)/6;
    
    p_i_f=p_i_f+(fi0_pi+4*fi1_pi+fi2_pi)/6;
    
    d_r_f=d_r_f+(fi0_dr+4*fi1_dr+fi2_dr)/6;
    
    e_r=ifft(e_r_f); e_i=ifft(e_i_f); p_r=ifft(p_r_f); p_i=ifft(p_i_f); d_r=ifft(d_r_f);
    
    norm_vari=norm([e_r; e_i; p_r; p_i; d_r],'inf');
    
    PP = log ( norm_vari );
    
    SP = SP + PP;
    
    e_r = e_r / norm_vari; e_i = e_i / norm_vari; p_r = p_r / norm_vari; p_i = p_i / norm_vari; d_r = d_r / norm_vari;

    if (mod(m,1e5)==0)
       disp(m);
       disp(m * tau);
       disp( SP / ( m * tau ) );
    end   

end

SP2 = SP / ( Nt_end * tau );

toc;











