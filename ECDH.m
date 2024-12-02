function [k2_x,k2_y,k1_x,k1_y]=ECDH(x,y)
a=4;
b=20;
p=997;

G_x = 40;
G_y = 37;

n_alice = x;     % alice's private key
[p_alice_x,p_alice_y]=add(G_x,G_y,G_x,G_y,a,b,p);
for i = 1:n_alice-2
    [p_alice_x,p_alice_y]=add(p_alice_x,p_alice_y,G_x,G_y,a,b,p);
end

n_bob = y;       % bob's private key
[p_bob_x,p_bob_y]=add(G_x,G_y,G_x,G_y,a,b,p);
for i = 1:n_bob-2
    [p_bob_x,p_bob_y]=add(p_bob_x,p_bob_y,G_x,G_y,a,b,p);
end 

[k1_x,k1_y]=add(p_bob_x,p_bob_y,p_bob_x,p_bob_y,a,b,p);
for i = 1:n_alice-2
   [k1_x,k1_y] = add(k1_x,k1_y,p_bob_x,p_bob_y,a,b,p);
end

[k2_x,k2_y]=add(p_alice_x,p_alice_y,p_alice_x,p_alice_y,a,b,p);
for i = 1:n_bob-2
   [k2_x,k2_y] = add(k2_x,k2_y,p_alice_x,p_alice_y,a,b,p);
end

function [x3,y3]=add(x1,y1,x2,y2,a,b,p)
equalma=panduanequal(x1,y1,x2,y2);
if equalma==0
   lamdafenzi=mod(y2-y1,p);
   lamdafenmu=mod(x2-x1,p);
else
   lamdafenzi=mod(3*x1^2+a,p);
   lamdafenmu=mod(2*y1,p);
end
lamdafenmuniyuan=exgcd(p,lamdafenmu);
lamda=mod(lamdafenmuniyuan*lamdafenzi,p);
x3=lamda^2-x1-x2;
x3=mod(x3,p);
y3=lamda*(x1-x3)-y1;
y3=mod(y3,p);
end

%判断P是否为±Q
function equalma=panduanequal(x1,y1,x2,y2)
if y1~=y2&&y1~=-y2||x1~=x2
    equalma=0;%0代表P≠±Q
else
    equalma=1;%1代表P=±Q
end
end
%扩展欧几里得算法求逆元
function niyuan = exgcd(a,b)
r1=a;
r2=b;
s1=1;
s2=0;
t1=0;
t2=1;
while r2>0
    q=floor(r1/r2);
    r=r1-q*r2;
    r1=r2;
    r2=r;
    s=s1-q*s2;
    s1=s2;
    s2=s;
    t=t1-q*t2;
    t1=t2;
    t2=t;
end
if t1<0
    t1=mod(t1,a);
end
niyuan=t1;
end

end