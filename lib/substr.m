function s2 = substr(s1, sub1, sub2)
% return the substring between the last sub1 in s1 and the first sub2 in
% s1.

t1 = strfind(s1, sub1);
if isempty(t1)
    t1 = 1;
end
t2 = strfind(s1, sub2);
if isempty(t2)
    t2 = length(s1);
end
s2 = s1(t1(end) + 1 : t2);