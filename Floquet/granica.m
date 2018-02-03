function [del1 Ne1]=granica(delta,N) %%%%для односвязных областей
% X=1;
% Q=[];
% for q=2:numel(x) %%%% всю первую строчку записываем как крайнюю в предположении что х всегда не убывает, а у в пределах одинаковых х всегда растет
%     if x(q)==x(1), X=[X q];else q1=q+1;break;end;
% end
% for q=q1:numel(x)-1 %%%%% только в предположении что у постоянно растет в пределах одинаковых х
%     if x(q)==x(q+1), Q=[Q q];elseif x(q+1)<x(numel(x)), X=[X q q+1];Q=x(q+1);else X=[X q+1 numel(x)];end
% end


p=0;
b=numel(delta);
for i=1:b
    if delta(i)==delta(1)
        p=p+1;delta1(p)=delta(i);N1(p)=N(i);
    elseif delta(i)<delta(i-1) %&& N(i-1)<300
        p=p+1;delta1(p)=delta(i-1);N1(p)=N(i-1);
    end;
end;
for i=1:p
    del1(i)=delta1(p+1-i);Ne1(i)=N1(p+1-i);
end;
for i=2:b
    if delta(i)<delta(i-1)
        del1(size(del1,2)+1)=delta(i);Ne1(size(Ne1,2)+1)=N(i);
    end;
end;
