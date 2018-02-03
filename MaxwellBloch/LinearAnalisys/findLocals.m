function x = findLocals(in)

side = 5;
isEdgeAvailable=0;

difs = diff(in);
x=[];
for i=2:numel(difs)
   if difs(i-1)*difs(i) < 0 && difs(i-1) > 0
       x = [x i];
   end
end