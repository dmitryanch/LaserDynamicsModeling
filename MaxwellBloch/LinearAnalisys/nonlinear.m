clear
clc
global sigma gamma a q k
%%%%%%%%%%%%%%%% A
% sigma=0.01;
% gamma=0.2;
% a=0.01;
% r=30;
% delta=-3;
% K=0:0.1:20;    
% lim=40;dq=0.01;
% nth=1;nimag=2;
%%%%%%%%%%%%% B
sigma=0.05;
gamma=0.0005;
a=0.01;
r=10;
delta=-1;
K=0:0.01:5;    
lim=10;dq=0.01;
nth=2;nimag=0.2;
%%%%%%%%%%%%%% C
% sigma=2.2;
% gamma=1;
% a=0.01;
% r=50;
% delta=-0.5;
% K=0:0.1:40;    
% lim=50;dq=0.1;
% nth=2;nimag=0.5;
% sigma=5;
% gamma=0.1;
% a=0.01;
% r=30;
% delta=-2;
% K=0:0.1:40;    
% lim=50;dq=0.1;
% nth=2;nimag=1;

gcaFontSize=18;
labelSize=22;
filename=['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)];

qq=-lim:dq:lim;
q=(-lim-max(K)):dq:(lim);
dq=abs(q(2)-q(1));
figure;
hold on;grid on;
x=[];y=[];
IMAG=zeros(numel(qq),numel(k));
%%
for i=1:numel(K)
    k=K(i);
    [x1 y1]=EigProgDelNegRandK(delta,r);
    
    %%% first
    locals=findLocals(x1);
    locals=locals(x1(locals)>0);
    plot(k,k+q(locals),'.','color',[0.4 0.4 0.4],'LineWidth',2)
    
    %%% second
%     locals=locals(imag(y1(locals))>1e-10);
%     ynorm=x1(locals)./((1+((a*(k+q(locals)).^2-delta)/(1+sigma)).^2).^1);
%     locals=locals(x1(locals)./((1+((a*(k+q(locals)).^2-delta)/(1+sigma)).^1).^2)>max(ynorm)*0.5);
%     if k<2.45 up=-1; else up=1e-1;    end
%     locals=locals(q(locals)+k<=up);  %% dirty hack - dont do like this
    [m ind]=max(x1(locals)./((1+((a*(k+q(locals)).^2-delta)/(1+sigma)).^2).^nth).*(imag(y1(locals))/max(imag(y1(locals)))*max(real(y1(locals)))).^nimag);x=[x k];y=[y q(locals(ind))+k];
%     
    %%% third
    tmp=x1./((1+((a*(k+q).^2-delta)/(1+sigma)).^2).^nth);
%     tmp=x1./((1+((a*(k+q).^2-delta)/(1+sigma)).^2).^nth).*(imag(y1)/max(imag(y1))*max(real(y1))).^nimag;
    start=find(round(q/dq)==round((-lim-k)/dq));finish=start+numel(qq)-1;
    IMAG(:,i)=real(tmp(start:finish).^0.5);

    drawnow;
end
%%
plot(K,K,'--k',K,-K,'--k',x,y,'.k','markers',15)
set(gca,'FontSize',gcaFontSize);
xlabel('k','FontSize',labelSize,'FontWeight','bold');ylabel('q','FontSize',labelSize,'FontWeight','bold');
    saveas(gcf,['locals ',filename,'.jpg'],'jpg');%close gcf;
%%
% figure;plot(x,y,'.black');
% hold on;grid on;plot(K,K,'--k',K,-K,'--k');xlim([min(x) max(x)]);ylim([-max(K) max(y)*1.2]);
% set(gca,'FontSize',gcaFontSize);
% xlabel('k','FontSize',18,'FontWeight','bold');ylabel('q','FontSize',18,'FontWeight','bold');
% % ylim([-50 20])    
% saveas(gcf,['maximums ',filename,'.jpg'],'jpg');%close gcf;
%%
figure;imagesc(K,qq,(((IMAG).^2)));
set(gca,'FontSize',gcaFontSize);
xlabel('k','FontSize',labelSize,'FontWeight','bold');ylabel('q','FontSize',labelSize,'FontWeight','bold');
colorbar('vert');colormap(gray)
saveas(gcf,['IMAG ',filename,'.jpg'],'jpg');%close gcf;

%%
% k=10;
% [x1 y1]=EigProgDelNegRandK(delta,r);
% % nth=1.2;nimag=1.9;
% locals=findLocals(x1);
% locals=locals(x1(locals)>0);
% q(locals)+k
% figure;plot(q+k,x1,q+k,imag(y1)/max(imag(y1))*max(real(y1))*0.5)
% figure;plot(q+k,y1(:)./((1+((a*(k+q(:)).^2-delta)/(1+sigma)).^1).^nth))
% figure;plot(q+k,(x1./((1+((a*(k+q).^2-delta)/(1+sigma)).^2).^nth).*(imag(y1)/max(imag(y1))*max(real(y1))).^nimag).^0.5)