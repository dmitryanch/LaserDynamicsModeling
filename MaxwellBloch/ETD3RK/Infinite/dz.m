function Z=dz(E)
Z=ft(ft(E).').';
% Z=ft(E);
Z=abs(Z).^2; %для эл. поля
% Z=abs(Z); %для инверсии населенностей