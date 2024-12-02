% n 分子  d 分母   m 模数
function y = modfrac( n,d,m )

n=mod(n,m);
d=mod(d,m);


i=1;
while mod(d*i,m) ~=1
    i=i+1;
end

y=mod(n*i,m);
end