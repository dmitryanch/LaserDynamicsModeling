% clear
% clc

% X=load('C:\Users\DEXP\Documents\Visual Studio 2013\Projects\CUDA_Diagram\floke s2 g0.1 nontrivial 200x201.txt');
X=load('C:\Users\DEXP\Documents\Visual Studio 2013\Projects\ModulationDiagram\ModulationDiagram\many 1d fixed M sin sigma=0.005 gamma=0.0005 a=0.0001 r=10 delta=-0.1 dt=0.5 Nt=2e6 Ndel=2e6 nt.txt');
% X=load('100x101 delta=-0.25 extT.txt');
minX=min(X(:,3)),maxX=max(X(:,3))
% minX=min(X),maxX=max(X)
aveX=(minX+maxX)/2

minK=min(X(:,4))
maxK=max(X(:,4))

%%
level=0.3;
lim=aveX+(maxX-aveX)*level
x=[];y=[];x1=[];y1=[];
separator1=1.4;
separator2=1.8;
for i=1:size(X,1),
    if X(i,3)>lim,
        if X(i,2) < separator1
            x=[x X(i,1)];y=[y X(i,2)];
        elseif X(i,2) > separator2
            x1=[x1 X(i,1)];y1=[y1 X(i,2)];
        end
    end;
end
[bx1 by2]=limits((x),(y));figure;
patch(by2,bx1,[0.9 0.9 0.9],'LineWidth',1);
hold on
[bx3 by4]=limits(x1,y1);
patch(by4,bx3,[0.9 0.9 0.9],'LineWidth',1);
xlim([0 3])
ylim([0 1])
filename=['100x100 -0.2 '];
set(gca,'FontSize',18);grid on;
xlabel('\Omega/\omega_r_e_l','FontSize',22,'FontWeight','bold');ylabel('m','FontSize',18,'FontWeight','bold');
% saveas(gcf,['diagr ',filename,'.jpg'],'jpg');%close gcf;
%% floke
level=0.46;
lim=aveX+(maxX-aveX)*level
x=[];y=[];x1=[];y1=[];
for i=1:size(X,1),
    if X(i,3)>lim,
        x=[x X(i,1)];y=[y X(i,2)];
    end;
end

[bx3 by4]=limits(x,y);
patch(by4,bx3,[0.4 0.4 0.4],'LineWidth',1);
xlim([-5 0])
ylim([0 50])
filename=['200x200 -0.25 svkrestin'];
set(gca,'FontSize',16);
xlabel('\delta','FontSize',18,'FontWeight','bold');
ylabel('r','FontSize',18,'FontWeight','bold');
% saveas(gcf,['qqqlin+diagr ',filename,'.jpg'],'jpg');%close gcf;
%%
N1=100;N2=5;
IMAG=-1*ones(N1,N2);
for i=1:N1
    for j=1:N2
        IMAG(i,j)=X((i-1)*N2+j,3);
%         IMAG(i,j)=X((i-1)*N1+j);
    end
end
figure;imagesc(flip(X(1:N2,2)),X(1:N2:size(X,1),1),flip(IMAG,2))
% figure;imagesc(W,M,IMAG)
xlabel('\Omega/\omega_r_e_l','FontSize',22,'FontWeight','bold');
ylabel('m','FontSize',18,'FontWeight','bold');
colorbar('vert')
filename=[num2str(N1),'x',num2str(N2),' cuda'];
% saveas(gcf,['imag ',filename,'.jpg'],'jpg');%close gcf;

%% by 1d fixed W
N1=100;N2=5;
for i=1:N1
    for j=1:N2
        IMAG1d(i,j)=X((i-1)*N2+j,3);
    end
end

for i=1:N2
    SP=IMAG1d(:,i);
    M=X(1:N2:size(X,1),1);
    figure;plot(M,SP) 
    xlabel('m','FontSize',22,'FontWeight','bold');
    ylabel('\mu','FontSize',18,'FontWeight','bold');
    title(['\Omega/\omega_r_e_l=',num2str(X(i,2))])
    grid on
    hold on 
    plot(M, zeros(size(M)),'k')
    saveas(gcf,['imag 1d sin nt w=',num2str(X(i,2)),'.jpg'],'jpg');close gcf;
end

%% by 1d fixed M
N1=5;N2=100;
for i=1:N1
    for j=1:N2
        IMAG1d(i,j)=X((i-1)*N2+j,3);
    end
end

for i=1:N1
    SP=IMAG1d(i,:);
    W=X(1:size(X,1)/N1,2);
    figure;plot(W,SP) 
    xlabel('\Omega/\omega_r_e_l','FontSize',22,'FontWeight','bold');
    ylabel('\mu','FontSize',18,'FontWeight','bold');
    title(['m=',num2str(X(i*N2,1))])
    grid on
    hold on 
    plot(W, zeros(size(W)),'k')
    saveas(gcf,['imag 1d fixed M sin nt m=',num2str(X(i*N2,1)),'.jpg'],'jpg');close gcf;
end

%% Dual figure
X1=load('C:\Users\DEXP\Documents\Visual Studio 2013\Projects\ModulationDiagram\ModulationDiagram\many 1d fixed M sin sigma=0.005 gamma=0.0005 a=0.0001 r=10 delta=-0.1 dt=0.5 Nt=2e6 Ndel=2e6 nt.txt');
X2=load('C:\Users\DEXP\Documents\Visual Studio 2013\Projects\ModulationDiagram\ModulationDiagram\many 1d fixed M cos sigma=0.005 gamma=0.0005 a=0.0001 r=10 delta=-0.1 dt=0.5 Nt=2e6 Ndel=2e6 nt.txt');
N1=5;N2=100;
for i=1:N1
    for j=1:N2
        IMAG1dsin(i,j)=X1((i-1)*N2+j,3);
        IMAG1dcos(i,j)=X2((i-1)*N2+j,3);
    end
end

for i=1:N1
    SPsin=IMAG1dsin(i,:);
    SPcos=IMAG1dcos(i,:);
    W=X1(1:size(X1,1)/N1,2);
    figure;plot(W,SPsin,'red');hold on;o=plot(W,SPcos,'blue') 
    xlabel('\Omega/\omega_r_e_l','FontSize',22,'FontWeight','bold');
    ylabel('\mu','FontSize',18,'FontWeight','bold');
    title(['m=',num2str(X1(i*N2,1))])
    grid on
    plot(W, zeros(size(W)),'k')
    legend('sin','cos','Location','southeast')
    saveas(gcf,['imag 1d fixed M sin+cos nt m=',num2str(X1(i*N2,1)),'.jpg'],'jpg');close gcf;
end