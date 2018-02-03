clear
clc
global sigma gamma a q k
sigma=0.1;
gamma=0.01;
a=0.01;
q=-50:0.01:50;
k=0;
r=20;
delta=-2;
filename=['sigma=',num2str(sigma),' gamma=',num2str(gamma),' a=',num2str(a),' r=',num2str(r),' delta=',num2str(delta)];
[x1 y1]=EigProgDelNegRandK(delta,r);
plot(q,real(y1),'black','LineWidth',2)
hold on
plot(q,imag(y1)/100,'--black','LineWidth',1)
xlim([-10 10])
plot(q,zeros(size(q)),'black','LineWidth',.025)
xlabel('q','FontSize',14,'FontWeight','bold')
grid on
saveas(gcf,['q1 ',filename,'.jpg'],'jpg');
[m ind1]=max(real(y1));
dk1=q(ind1)

%%
k=dk1
[x2 y2]=EigProgDelNegRandK(delta,r);
plot(q,real(y2),'black','LineWidth',2)
hold on
plot(q,imag(y2)/10,'--black','LineWidth',1)
xlim([-10 10])
plot(q,zeros(size(q)),'black','LineWidth',.025)
xlabel('q','FontSize',14,'FontWeight','bold')
grid on
% saveas(gcf,['q1 ',filename,'.jpg'],'jpg');
[m ind2]=max(real(y2));
dk2=q(ind2)