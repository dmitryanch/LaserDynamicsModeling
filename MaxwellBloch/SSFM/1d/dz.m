function Z=dz(E)
Z=ft(ft(E).').';
% Z=ft(E);
Z=abs(Z).^2; %��� ��. ����
% Z=abs(Z); %��� �������� �������������