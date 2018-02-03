function Y=f3(E,P,D)
global gamma r
Y=-gamma*(D-r+0.5*(conj(E)*P+E*conj(P)));