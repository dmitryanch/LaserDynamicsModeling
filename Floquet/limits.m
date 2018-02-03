function [x, y]=limits(X, Y)
% conditions: X is non decreased, Y is non decreased for same X values
x=X(1);y=Y(1);
x1=[];y1=[];
for i=2:numel(X)
    if X(i)>X(i-1)
       x1=[x1 X(i-1)];
       y1=[y1 Y(i-1)];
       x=[x X(i)];
       y=[y Y(i)];
    elseif i==numel(X)
        x1=[x1 X(i)];
       y1=[y1 Y(i)];
    end
end
x=[x flip(x1)];
y=[y flip(y1)];