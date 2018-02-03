clear
clc
global sigma gamma a q k
sigma=2;
gamma=0.1;
a=0.01;
lim=40;
q=-lim:0.02:lim;
r=10;
delta=-0.5;
filename=['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)];
    
for k=0%[1,3,5,9,11,13,15]
    [x1 y1]=EigProgDelNegRandK(delta,r);
%     y1=y1./((1+((a*(k+q).^2-delta)/(1+sigma)).^2));
%     y1=x1./((1+((a*(k+q).^2-delta)/(1+sigma)).^2)).*(imag(y1)/max(imag(y1))*max(real(y1))).^2;[m ind]=max(y1);q(ind)+k
    figure;plot(q,real(y1),'black','LineWidth',2)
    hold on
    plot(q,imag(y1)/max(imag(y1))*max(real(y1))*0.5,'--k','LineWidth',1)
%     xlim([-10 10])
%     ylim([-0.1 0.1])
    plot(q,zeros(size(q)),'k','LineWidth',.025)
    xlabel('q','FontSize',18,'FontWeight','bold')
    grid on
%     saveas(gcf,['k=',num2str(k),' ',filename,'.jpg'],'jpg');%close gcf;
end
