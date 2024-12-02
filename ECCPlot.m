%a,b分别为椭圆的参数，p为一个质数
function [x,y] =  ECCPlot( a,b,p)

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
plot(x,y,'*')
hold on
grid on
end
