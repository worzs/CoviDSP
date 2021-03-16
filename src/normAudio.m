function s_norm = normAudio(s)
%{
 In:
 s - Audio vector 

 Out:
 s_norm - audio vector, with values between -1 and 1. 
%}

s_norm = (s-mean(s))/ max(abs(s));

cond = abs(s_norm) > db2mag(-10); 
index_begin = find(cond, 1, 'first');
index_end = find(cond, 1, 'last');
s_norm = s_norm(index_begin:index_end);

end
