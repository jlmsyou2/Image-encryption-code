%完整代码将在论文接收后上传
function [c2x,c2y] = ECC(key)
a=4;
b=20;
p=997;
% disp(['您输入的a,b,p依次为：',num2str(a),',',num2str(b),',',num2str(p)]);
% if mod(4*a^2+27*b^3,p)~=0
%    disp(['您输入的是F',num2str(p),'上的椭圆曲线：y^2=x^3+',num2str(a),'x+',num2str(b)]);
% else
%    disp('您输入的参数有误，不能构成椭圆曲线！请重新输入');
%    a=input('请输入一个a值：');
%    b=input('请输入一个b值：');
%    p=input('请输入一个p值：');
% end
% [x,y]=evaldian(a,b,p);
x0 = 40;
y0 = 37;
n=jieshu(x0,y0,a,b,p);
[xq,yq]=add(x0,y0,x0,y0,a,b,p);
a0 = 10;
for i1=1:a0-2
    [xq,yq]=add(x0,y0,xq,yq,a,b,p);
end
disp(['计算得到的公钥Q为：(',num2str(xq),',',num2str(yq),')']);
m=key;
R=20;
[c1x,c1y,c2x,c2y]=jiami(m,R,x0,y0,xq,yq,a,b,p);
disp(['加密后的密文(c1,c2)=(',num2str(c1x),',',num2str(c1y),',',num2str(c2x),',',num2str(c2y),')']);
jiemim=jiemi(a0,c1x,c1y,c2x,a,b,p);
disp(['解密得到的密文m=',num2str(jiemim)]);


function [x,y]=evaldian(a,b,p)
%计算出0-p范围内所有符合椭圆曲线的点
x=[];
y=[];
index = 1;
for xr = 0:1:p
    mm = mod(xr^3+a*xr+b ,p);    
    for yr=0:1:p
        if mod(yr^2,p) == mm
            x(index)=xr;
            y(index)=yr;
            index = index+1;
        end
    end
end
end

%加法及倍点运算
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

%计算G的阶数
function n=jieshu(x0,y0,a,b,p)
cnt=2;
x0fu=x0;
y0fu=mod(-y0,p);
disp(['(',num2str(x0),',',num2str(y0),')']);
[x1,y1]=add(x0,y0,x0,y0,a,b,p);
disp(['(',num2str(x1),',',num2str(y1),')']);
while x1~=x0fu||y1~=y0fu
    cnt=cnt+1;
    [x1,y1]=add(x1,y1,x0,y0,a,b,p);
    disp(['(',num2str(x1),',',num2str(y1),')']);
end
n=cnt+1;
end

%加密算法
function [c1x,c1y,c2_x,c2_y]=jiami(m,r,xp,yp,xq,yq,a,b,p)
[c1x,c1y]=add(xp,yp,xp,yp,a,b,p);
for i=1:r-2
    [c1x,c1y]=add(c1x,c1y,xp,yp,a,b,p);
end
[rqx,rqy]=add(xq,yq,xq,yq,a,b,p);
for i=1:r-2
    [rqx,rqy]=add(rqx,rqy,xq,yq,a,b,p);
end

c2_x=mod(m*rqx,p);
c2_y=mod(m*rqy,p);

% [c2_x,c2_y]=add(rqx,rqy,rqx,rqy,a,b,p);
% for i=1:m-2
%     [c2_x,c2_y]=add(c2_x,c2_y,rqx,rqy,a,b,p);
% end
end

%解密算法
function jiemim=jiemi(a0,c1x,c1y,c2,a,b,p)
[x2,y2]=add(c1x,c1y,c1x,c1y,a,b,p);
for i=1:a0-2
    [x2,y2]=add(x2,y2,c1x,c1y,a,b,p);
end
x2niyuan=exgcd(p,x2);
jiemim=mod(c2*x2niyuan,p);
end
end
