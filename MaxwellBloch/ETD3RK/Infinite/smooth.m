function y=smooth(x,k)
[m n]=size(x);
s=3; %по квадрату со стороной s*2+1 точек
x=[x(m-s+1:m,:);x;x(1:s,:)];
x=[x(:,n-s+1:n),x,x(:,1:s)];
y=x;
for l=1:k
for i=1+s:size(x,1)-s
    for j=1+s:size(x,2)-s
        y(i,j)=sum(sum(x(i-s:i+s,j-s:j+s)))/(2*s+1)^2;
    end
end
y=y(s+1:size(y,1)-s,s+1:size(y,2)-s);
x=y;
x=[x(m-s+1:m,:);x;x(1:s,:)];
x=[x(:,n-s+1:n),x,x(:,1:s)];
y=x;
end
y=y(s+1:size(y,1)-s,s+1:size(y,2)-s);