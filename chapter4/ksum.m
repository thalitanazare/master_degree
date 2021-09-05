function s = ksum(x)
[~, idx] = sort(abs(x),'descend'); 
x = x(idx);
s=0; e=0;
for i=1:numel(x)
   s_old = s; 
   s = s + x(i); 
   e = e + x(i) - (s - s_old);  
end
s = s + e;
return