function Y=ft(y)
global x0 x dx
m=size(y,1);
%%%%dx=1/m;
% if m==1, x0=0; else x0=linspace(-d/2,d/2,m); end; %���������� � ������� ����
% x0=(0:m-1)*dx;
% x=-xdz*x0; %���������� � ������� ����
% % % dx=abs(x0(1)-x0(2));
% x=xdz*linspace(-max(x0)/2,max(x0)/2,m);
lambda=1e-2;%����� �����
f=1e2;% �������� ����������
[X0 X]=meshgrid(x0,x);
S=exp(2*pi*1i*X0.*X/lambda/f);
S1=1i*exp(-1i*4*pi*f/lambda)/f/lambda;
for m=1:size(y,2)
    y1=meshgrid(y(:,m));
    Y(:,m)=sum(dx*y1.*S,2)*S1;
end

% ���������� �����-��������������
% % % [X0 X]=meshgrid(1:size(y,1));
% % % for m=1:size(y,2)
% % %     y1=meshgrid(y(:,m));
% % %     Y(:,m)=sum(y1.*exp(-2*pi*i*(X0-1).*(X-1)/size(y,1)),2);
% % % end