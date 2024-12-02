clear; 
clc;
%-------------------------------------------------------------------------%
%------------------------读取图像并通过hash函数产生对应哈希值---------------%
% img=imread('medical_chest_X_ray_png2.png');
% img=imread('medical_CT_COVID_png3.png');
img=imread('medical_HAM_jpg1.jpg');
% img=imread('ISIC_0024306.jpg');
I_r=uint8(img(:,:,1));
I_g=uint8(img(:,:,2));
I_b=uint8(img(:,:,3));
[M,N]=size(I_r);  
 
I_avg_hash=DataHash(img,'array','SHA-256');
I_avg_hash_value = reshape(I_avg_hash,[64,1]);%转换为1*64矩阵
I_avg_hash_dec = hex2dec(I_avg_hash_value);%转为dec（十进制）

a0 = sum(I_avg_hash_dec(1:12));
b0 = sum(I_avg_hash_dec(13:24));
x0 = sum(I_avg_hash_dec(25:36));
y0 = sum(I_avg_hash_dec(37:48));
Sr = sum(I_avg_hash_dec(49:56));
Sg = sum(I_avg_hash_dec(53:60));
Sb = sum(I_avg_hash_dec(57:64));

[a1_x,a1_y]=ECC(a0);
[b1_x,b1_y]=ECC(b0);
[x1_x,x1_y]=ECC(x0);
[y1_x,y1_y]=ECC(y0);
[Sr1_x,Sr1_y]=ECC(Sr);
[Sg1_x,Sg1_y]=ECC(Sg);
[Sb1_x,Sb1_y]=ECC(Sb);
[Skey1_x,Skey1_y,Skey2_x,Skey2_y]=ECDH(71,59);

a2 = mod(log2(a0*Skey1_x)+(a1_x*b1_y),M);
b2 = mod(log2(b0*Skey1_y)+(b1_x*a1_y),N);
x2 = sin(pi*((x0*Skey1_x)/x1_x+y1_x));
y2 = sin(pi*((y0*Skey1_y)/x1_y+y1_y));
Sr2= mod(abs(Sr1_x-Sg1_y),M);
Sg2= mod(abs(Sg1_x-Sb1_y),M);
Sb2= mod(abs(Sb1_x-Sr1_y),M);

a = Sr2;
b = Sg2;
c = Sb2;
A = a+b+c;
%-------------------------------------------------------------------------%
%------------------------迭代出加密所需的混沌序列---------------------------%
x=x2;
y=y2;
for i=1:1000            %step1 去除前1000次迭代结果
      x = mod(10*pi*((1-y)*sin(a2*sin(pi*x)*b2*y*(1-y^2))+b2*y^2),1);
      y = mod(10*pi*((1-x)*sin(a2*sin(pi*y)*b2*x*(1-x^2))+a2*x^2),1);
end

H1=zeros(1,A+256 );
H2=zeros(1,A+256 );
H1(1)=x;
H2(1)=y;
for r=1:A+256-1  %step2 迭代M*N-1次,得到的序列值被用于决定要抽取的像素的坐标、S_box及其横纵坐标规则   
    H1(r+1) = mod(10*pi*((1-H2(r))*sin(a2*sin(pi*H1(r))*b2*H2(r)*(1-H2(r)^2))+b2*H2(r)^2),1);
    H2(r+1) = mod(10*pi*((1-H1(r))*sin(a2*sin(pi*H2(r))*b2*H1(r)*(1-H1(r)^2))+a2*H1(r)^2),1);
end
%-------------------构造被抽取出的像素值的坐标和S_box_A---------------------%
p1 = mod(floor(H1(1:A).*10^13),M)+1; %被抽取出的像素值横坐标
p2 = mod(floor(H2(1:A).*10^13),N)+1; %被抽取出的像素值纵坐标
S_box_H1 = H1(A+1:A+256);
[~,S_box_idx] = sort(S_box_H1);  
S_box_A = reshape(S_box_idx,[16,16])-1;%S_box0~255替换值
%-------------------------------------------------------------------------%
H3=zeros(1,M*N);
H4=zeros(1,M*N);
H3(1)=H1(A+256);
H4(1)=H2(A+256);
for r=1:M*N-1           %step2 迭代M*N-1次，得到的序列被用于获得S_box_B坐标规则
    H3(r+1) = mod(10*pi*((1-H4(r))*sin(a2*sin(pi*H3(r))*b2*H4(r)*(1-H4(r)^2))+b2*H4(r)^2),1);
    H4(r+1) = mod(10*pi*((1-H3(r))*sin(a2*sin(pi*H4(r))*b2*H3(r)*(1-H3(r)^2))+a2*H3(r)^2),1);
end
%-------------------------------------------------------------------------%
%------------------------下述步骤将每个平面被取出的像素依次放入S序列中-------%
S_r=zeros(1,a);
S_g=zeros(1,b);
S_b=zeros(1,c);
S=zeros(1,A); 
for r=1:a
    S_r(r)=I_r(p1(r),p2(r));
    S(r)=S_r(r);
end
for r=a+1:a+b
    S_g(r-a)=I_g(p1(r),p2(r));
    S(r)=S_g(r-a);
end
for r=a+b+1:a+b+c
    S_b(r-a-b)=I_b(p1(r),p2(r));
    S(r)=S_b(r-a-b);
end
S_idx = p1.*p2;
[~, idx] = sort(S_idx);% 使用S_sort的排序顺序对S进行排序,以此打乱每个平面取出的像素，然后得到S_idx。
S_sort = S(idx);
%-------------------------------------------------------------------------%
%------------------------对抽取出来的A=a+b+c个像素值进行S_box替换-----------%
S_sort_H = mod(S_sort,16)+1;
S_sort_L = floor(S_sort/16)+1; 
S_sort_box = zeros(1,A);
          for i=1:A          
              S_sort_box(i)=S_box_A(S_sort_H(i),S_sort_L(i));
          end
%-------------------------------------------------------------------------%
%-------------------获取被取出像素后的三个平面I1_r,I1_g,I1_b----------------%
I1_r=zeros(1,M*N-a);
I1_g=zeros(1,M*N-b);
I1_b=zeros(1,M*N-c);

D=p1+(p2-1)*M; %D和每个平面被取出来的位置一一对应，长度为a+b+c。
Da=D(1:a);%Da和r平面被取出像素坐标对应。
Db=D(a+1:a+b);%Db和g平面被取出像素坐标对应。
Dc=D(a+b+1:a+b+c);%Dc和b平面被取出像素坐标对应。

D1a=sort(Da);%对Da排序，方便从I_r中按出现顺序提取出像素，下同
D1b=sort(Db);
D1c=sort(Dc);

for i=1:a+1 %去除I_r的中的a个像素
    if(i==1)
        I1_r(1:D1a(i)-1)=I_r(1:D1a(i)-1);
    elseif(i<=a)
        z=D1a(i-1);
        I1_r(z-i+2:D1a(i)-i)=I_r(z+1:D1a(i)-1);
    else
        z=D1a(i-1);
        I1_r(z-a+1:M*N-a)=I_r(z+1:M*N);
    end
end
for i=1:b+1 %去除I_g的中的b个像素
    if(i==1)
        I1_g(1:D1b(i)-1)=I_g(1:D1b(i)-1);
    elseif(i<=b)
        z=D1b(i-1);
        I1_g(z-i+2:D1b(i)-i)=I_g(z+1:D1b(i)-1);
    else
        z=D1b(i-1);
        I1_g(z-b+1:M*N-b)=I_g(z+1:M*N);
    end
end
for i=1:c+1 %去除I_b的中的c个像素
    if(i==1)
        I1_b(1:D1c(i)-1)=I_b(1:D1c(i)-1);
    elseif(i<=c)
        z=D1c(i-1);
        I1_b(z-i+2:D1c(i)-i)=I_b(z+1:D1c(i)-1);
    else
        z=D1c(i-1);
        I1_b(z-c+1:M*N-c)=I_b(z+1:M*N);
    end
end
%-------------------------------------------------------------------------%
%-------------------------------将I1_r,I1_g,I1_b重塑为[m,n]矩阵-----------%
s=M*N-a;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
lr = length(F);
mr = F(ceil(lr/2));
nr = s/mr;
I1_r = reshape(I1_r,[mr,nr]);

s=M*N-b;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
lg = length(F);
mg = F(ceil(lg/2));
ng = s/mg;
I1_g = reshape(I1_g,[mg,ng]);

s=M*N-c;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
lb = length(F);
mb = F(ceil(lb/2));
nb = s/mb;
I1_b = reshape(I1_b,[mb,nb]);
%-------------------------------------------------------------------------%
%--------------------------------设置S_box_B参数值,6种坐标规则--------------%

G1=DataHash(S_sort_box(1:a),'array','SHA-256');
G1_hash_value = reshape(G1,[32,2]);%转换为1*64矩阵
G1_hash_dec = hex2dec(G1_hash_value);%转为dec（十进制）

G2=DataHash(S_sort_box(a+1:a+b),'array','SHA-256');
G2_hash_value = reshape(G2,[32,2]);%转换为1*64矩阵
G2_hash_dec = hex2dec(G2_hash_value);%转为dec（十进制）

G3=DataHash(S_sort_box(a+b+1:a+b+c),'array','SHA-256');
G3_hash_value = reshape(G3,[32,2]);%转换为1*64矩阵
G3_hash_dec = hex2dec(G3_hash_value);%转为dec（十进制）

S_box_H2 = H2(A+1:A+256);
[~,S_box_idx] = sort(S_box_H2);  
S_box_B = reshape(S_box_idx,[16,16])-1;%S_box0~255替换值

G1_rows = mod(G1_hash_dec,16);
G1_columns = floor(G1_hash_dec/16);
G2_rows = mod(G2_hash_dec,16);
G2_columns = floor(G2_hash_dec/16);
G3_rows = mod(G3_hash_dec,16);
G3_columns = floor(G3_hash_dec/16);

row1 = G1_rows(1:16);  row2 = G1_rows(17:32);
row3 = G2_rows(1:16);  row4 = G2_rows(17:32);
row5 = G3_rows(1:16);  row6 = G3_rows(17:32);
[~, rowidx1] = sort(row1);[~, rowidx2] = sort(row2);
[~, rowidx3] = sort(row3);[~, rowidx4] = sort(row4);
[~, rowidx5] = sort(row5);[~, rowidx6] = sort(row6);
[~, rowidx11] = sort(rowidx1);[~, rowidx22] = sort(rowidx2);
[~, rowidx33] = sort(rowidx3);[~, rowidx44] = sort(rowidx4);
[~, rowidx55] = sort(rowidx5);[~, rowidx66] = sort(rowidx6);

columns1 = G1_columns(1:16);  columns2 = G1_columns(17:32);
columns3 = G2_columns(1:16);  columns4 = G2_columns(17:32);
columns5 = G3_columns(1:16);  columns6 = G3_columns(17:32);
[~, columnsidx1] = sort(columns1);[~, columnsidx2] = sort(columns2);
[~, columnsidx3] = sort(columns3);[~, columnsidx4] = sort(columns4);
[~, columnsidx5] = sort(columns5);[~, columnsidx6] = sort(columns6);
[~, columnsidx11] = sort(columnsidx1);[~, columnsidx22] = sort(columnsidx2);
[~, columnsidx33] = sort(columnsidx3);[~, columnsidx44] = sort(columnsidx4);
[~, columnsidx55] = sort(columnsidx5);[~, columnsidx66] = sort(columnsidx6);

%-------------------------------------------------------------------------%
%-------------------------分解矩阵I1_r,I1_g,I1_b为上下两部分----------------%
I1_rL = floor(I1_r/16)+1;
I1_rH = mod(I1_r,16)+1;
I1_gL = floor(I1_g/16)+1;
I1_gH = mod(I1_g,16)+1;
I1_bL = floor(I1_b/16)+1;
I1_bH = mod(I1_b,16)+1;
%-------------------------------------------------------------------------%
%-------------------------------进行S_box_B替换得到I2_r,I2_g,I2_b----------%
pr5 = mod(floor(H3(1:M*N-a)*10^13),6)+1; %r平面选择横纵坐标的规则
pg5 = mod(floor(H4(1:M*N-b)*10^13),6)+1; %g平面选择横纵坐标的规则
pb5 = mod(floor(H3(1:M*N-c)*10^13),6)+1; %b平面选择横纵坐标的规则
I2_r = zeros(size(I1_rL));
I2_g = zeros(size(I1_gL));
I2_b = zeros(size(I1_bL));
for i=1:M*N-a
   switch(pr5(i))
       case 1
           I2_r(i) = S_box_B(rowidx11(I1_rH(i)),columnsidx11(I1_rL(i)));
       case 2
           I2_r(i) = S_box_B(rowidx22(I1_rH(i)),columnsidx22(I1_rL(i)));
       case 3
           I2_r(i) = S_box_B(rowidx33(I1_rH(i)),columnsidx33(I1_rL(i)));
       case 4
           I2_r(i) = S_box_B(rowidx44(I1_rH(i)),columnsidx44(I1_rL(i)));
       case 5
           I2_r(i) = S_box_B(rowidx55(I1_rH(i)),columnsidx55(I1_rL(i)));
       case 6
           I2_r(i) = S_box_B(rowidx66(I1_rH(i)),columnsidx66(I1_rL(i)));
   end
end
for i=1:M*N-b
   switch(pg5(i))
       case 1
           I2_g(i) = S_box_B(rowidx11(I1_gH(i)),columnsidx11(I1_gL(i)));
       case 2
           I2_g(i) = S_box_B(rowidx22(I1_gH(i)),columnsidx22(I1_gL(i)));
       case 3
           I2_g(i) = S_box_B(rowidx33(I1_gH(i)),columnsidx33(I1_gL(i)));
       case 4
           I2_g(i) = S_box_B(rowidx44(I1_gH(i)),columnsidx44(I1_gL(i)));
       case 5
           I2_g(i) = S_box_B(rowidx55(I1_gH(i)),columnsidx55(I1_gL(i)));
       case 6
           I2_g(i) = S_box_B(rowidx66(I1_gH(i)),columnsidx66(I1_gL(i)));
   end
end
for i=1:M*N-c
   switch(pb5(i))
       case 1
           I2_b(i) = S_box_B(rowidx11(I1_bH(i)),columnsidx11(I1_bL(i)));
       case 2
           I2_b(i) = S_box_B(rowidx22(I1_bH(i)),columnsidx22(I1_bL(i)));
       case 3
           I2_b(i) = S_box_B(rowidx33(I1_bH(i)),columnsidx33(I1_bL(i)));
       case 4
           I2_b(i) = S_box_B(rowidx44(I1_bH(i)),columnsidx44(I1_bL(i)));
       case 5
           I2_b(i) = S_box_B(rowidx55(I1_bH(i)),columnsidx55(I1_bL(i)));
       case 6
           I2_b(i) = S_box_B(rowidx66(I1_bH(i)),columnsidx66(I1_bL(i)));
   end
end
%-------------------------------------------------------------------------%
%-------------------------将I2_r,I2_g,I2_b重塑为1个矩阵I2_All--------------%
% I2_All = zeros(1,3*M*N-A);
I2_All(1:M*N-a)=I2_r;
I2_All(M*N-a+1:2*M*N-a-b)=I2_g;
I2_All(2*M*N-a-b+1:3*M*N-A)=I2_b;

s=3*M*N-A;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
l_I2_All = length(F);
m_I2_All = F(ceil(l_I2_All/2));
n_I2_All = s/m_I2_All;
I2_All = reshape(I2_All,[m_I2_All,n_I2_All]);
%-------------------------------------------------------------------------%
%--------------------------------------对矩阵I2_All进行扩散和旋转置乱-------%
%-------------------------------------------------------------------------%
H5=zeros(1,6*(M+N));
H6=zeros(1,6*(M+N));
H5(1)=H3(M*N);
H6(1)=H4(M*N);
for r=1:6*(M+N)-1           %step2 迭代M*N-1次，得到的序列被用于获得S_box_B坐标规则
    H5(r+1) = mod(10*pi*((1-H6(r))*sin(a2*sin(pi*H5(r))*b2*H6(r)*(1-H6(r)^2))+b2*H6(r)^2),1);
    H6(r+1) = mod(10*pi*((1-H5(r))*sin(a2*sin(pi*H6(r))*b2*H5(r)*(1-H5(r)^2))+a2*H5(r)^2),1);
end
Z1 = mod(ceil(H5(1:3*(M+N))*10^13),4)+1;
Z2 = mod(ceil(H5(3*(M+N)+1:6*(M+N))*10^13),A)+1;
hang = mod(ceil(H6(1:3*(M+N))*10^13),m_I2_All)+1;
lie = mod(ceil(H6(3*(M+N)+1:6*(M+N))*10^13),n_I2_All)+1;

H7=zeros(1,3*M*N-A);
H8=zeros(1,3*M*N-A);
H7(1)=H5(6*(M+N));
H8(1)=H6(6*(M+N));
for r=1:3*M*N-A-1           
    H7(r+1) = mod(10*pi*((1-H8(r))*sin(a2*sin(pi*H7(r))*b2*H8(r)*(1-H8(r)^2))+b2*H8(r)^2),1);
    H8(r+1) = mod(10*pi*((1-H7(r))*sin(a2*sin(pi*H8(r))*b2*H7(r)*(1-H7(r)^2))+a2*H7(r)^2),1);
end

I2_All_256 = mod(ceil(H7*10^13),256);
I2_All_1 = zeros(size(I2_All));
I2_All_1(1) = mod(I2_All(1)+I2_All_256(1),256); 
for i=2:3*M*N-A
    I2_All_1(i) = mod(I2_All_1(i-1)+I2_All(i)+I2_All_256(i),256); 
end
I2_All_2 = I2_All_1;
for i=1:3*(M+N)
    switch(Z1(i))
        case 1 
                I2_All_2(:, lie(i)) = circshift(I2_All_2(:, lie(i)), -Z2(i)); % 对H5(i)列进行向下循环移位Z7（j） 
        case 2 
                I2_All_2(hang(i), :) = circshift(I2_All_2(hang(i), :), -Z2(i)); % 对每H6(i)行进行左循环移位Z6（j）
        case 3 
                I2_All_2(:, lie(i)) = circshift(I2_All_2(:, lie(i)), Z2(i)); % 对H5(i)列进行向下循环移位Z7（j）
        case 4
                I2_All_2(hang(i), :) = circshift(I2_All_2(hang(i), :), Z2(i)); % 对每H6(i)行进行向右循环移位Z6（j）
    end 
end
%-------------------------------------------------------------------------%
%-----------------------------将矩阵I2_All重塑为I2_r,I2_g,I2_b-------------%
I2_r = I2_All_2(1:M*N-a);
I2_g = I2_All_2(M*N-a+1:2*M*N-a-b);
I2_b = I2_All_2(2*M*N-a-b+1:3*M*N-a-b-c);
I2_r = reshape(I2_r,[mr,nr]);
I2_g = reshape(I2_g,[mg,ng]);
I2_b = reshape(I2_b,[mb,nb]);
%-------------------------------------------------------------------------%
%---将a,b,c个元素分别放入I2_r,I2_g,I2_b中，得到[M,N]的矩阵I3_r,I3_g,I3_b----%
I3_r = zeros(M,N);
I3_g = zeros(M,N);
I3_b = zeros(M,N);
for i=1:a+1 %将S_idx中的a个元素插入到I2_r(Da(i))
    if(i==1)
        I3_r(1:D1a(i)-1)=I2_r(1:D1a(i)-1);
        I3_r(D1a(i))=S_sort_box(i);
    elseif(i<=a)
        I3_r(D1a(i-1)+1:D1a(i)-1)=I2_r(D1a(i-1)-i+2:D1a(i)-i);
        I3_r(D1a(i))=S_sort_box(i);
    else
        I3_r(D1a(a)+1:M*N)=I2_r(D1a(a)-a+1:mr*nr);
    end
end
for i=1:b+1 %将S_idx中的b个元素插入到I2_r(Da(i))
    if(i==1)
        I3_g(1:D1b(i)-1)=I2_g(1:D1b(i)-1);
        I3_g(D1b(i))=S_sort_box(a+i);
    elseif(i<=b)
        I3_g(D1b(i-1)+1:D1b(i)-1)=I2_g(D1b(i-1)-i+2:D1b(i)-i);
        I3_g(D1b(i))=S_sort_box(a+i);
    else
        I3_g(D1b(b)+1:M*N)=I2_g(D1b(b)-b+1:mg*ng);
    end
end
for i=1:c+1 %将S_idx中的c个元素插入到I2_r(Da(i))
    if(i==1)
        I3_b(1:D1c(i)-1)=I2_b(1:D1c(i)-1);
        I3_b(D1c(i))=S_sort_box(a+b+i);
    elseif(i<=c)
        I3_b(D1c(i-1)+1:D1c(i)-1)=I2_b(D1c(i-1)-i+2:D1c(i)-i);
        I3_b(D1c(i))=S_sort_box(a+b+i);
    else
        I3_b(D1c(c)+1:M*N)=I2_b(D1c(c)-c+1:mb*nb);
    end
end
I3(:,:,1)=uint8(I3_r);
I3(:,:,2)=uint8(I3_g);
I3(:,:,3)=uint8(I3_b);
figure(1)
imshow(I3);


%-----------------------------------解密-----------------------------------%
%---------从密文I3中抽取出a+b+c个元素，获取D3_r,D3_g,D3_b-------------------%

D3_r=zeros(1,M*N-a);
D3_g=zeros(1,M*N-b);
D3_b=zeros(1,M*N-c);

D=p1+(p2-1)*M; %D和每个平面被取出来的位置一一对应，长度为a+b+c。
Da=D(1:a);%Da和r平面被取出像素坐标对应。
Db=D(a+1:a+b);%Db和g平面被取出像素坐标对应。
Dc=D(a+b+1:a+b+c);%Dc和b平面被取出像素坐标对应。

D1a=sort(Da);%对Da排序，方便从I_r中按出现顺序提取出像素，下同
D1b=sort(Db);
D1c=sort(Dc);

for i=1:a+1 %去除D3_r的中的a个像素
    if(i==1)
        D3_r(1:D1a(i)-1)=I3_r(1:D1a(i)-1);
    elseif(i<=a)
        z=D1a(i-1);
        D3_r(z-i+2:D1a(i)-i)=I3_r(z+1:D1a(i)-1);
    else
        z=D1a(i-1);
        D3_r(z-a+1:M*N-a)=I3_r(z+1:M*N);
    end
end
for i=1:b+1 %去除D3_g的中的b个像素
    if(i==1)
        D3_g(1:D1b(i)-1)=I3_g(1:D1b(i)-1);
    elseif(i<=b)
        z=D1b(i-1);
        D3_g(z-i+2:D1b(i)-i)=I3_g(z+1:D1b(i)-1);
    else
        z=D1b(i-1);
        D3_g(z-b+1:M*N-b)=I3_g(z+1:M*N);
    end
end
for i=1:c+1 %去除D3_b的中的c个像素
    if(i==1)
        D3_b(1:D1c(i)-1)=I3_b(1:D1c(i)-1);
    elseif(i<=c)
        z=D1c(i-1);
        D3_b(z-i+2:D1c(i)-i)=I3_b(z+1:D1c(i)-1);
    else
        z=D1c(i-1);
        D3_b(z-c+1:M*N-c)=I3_b(z+1:M*N);
    end
end
%-------------------------------------------------------------------------%
%-------------------------------将D3_r,D3_g,D3_b重塑为[m,n]矩阵-----------%
s=M*N-a;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
lr = length(F);
mr = F(ceil(lr/2));
nr = s/mr;
D3_r = reshape(D3_r,[mr,nr]);

s=M*N-b;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
lg = length(F);
mg = F(ceil(lg/2));
ng = s/mg;
D3_g = reshape(D3_g,[mg,ng]);

s=M*N-c;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
lb = length(F);
mb = F(ceil(lb/2));
nb = s/mb;
D3_b = reshape(D3_b,[mb,nb]);

a1 = bitxor(I2_b,D3_b);
countA = sum(a1(:) ~= 0);

%-------------------------------------------------------------------------%
%-------------------------将D3_r,D3_g,D3_b重塑为1个矩阵D3_All--------------%
D3_All_2(1:M*N-a)=D3_r;
D3_All_2(M*N-a+1:2*M*N-a-b)=D3_g;
D3_All_2(2*M*N-a-b+1:3*M*N-A)=D3_b;

s=3*M*N-A;
F =[];
k=1;
for i=1:round(s/2) %找所有因子
	if (mod(s,i)==0)
		F(k)=i;
		k=k+1;
	end
end
l_D3_All_2 = length(F);
m_D3_All_2 = F(ceil(l_D3_All_2/2));
n_D3_All_2 = s/m_D3_All_2;
D3_All_2 = reshape(D3_All_2,[m_D3_All_2,n_D3_All_2]);
f1 = bitxor(I2_All_2,D3_All_2);
countF = sum(f1(:) ~= 0);
%-------------------------------------------------------------------------%
%-----------------------------------对D3_All进行旋转置乱逆过程--------------%
for i=3*(M+N):-1:1
    switch(Z1(i))
        case 1 
                D3_All_2(:, lie(i)) = circshift(D3_All_2(:, lie(i)), Z2(i)); % 对H5(i)列进行向下循环移位Z7（j） 
        case 2 
                D3_All_2(hang(i), :) = circshift(D3_All_2(hang(i), :), Z2(i)); % 对每H6(i)行进行左循环移位Z6（j）
        case 3 
                D3_All_2(:, lie(i)) = circshift(D3_All_2(:, lie(i)), -Z2(i)); % 对H5(i)列进行向下循环移位Z7（j）
        case 4
                D3_All_2(hang(i), :) = circshift(D3_All_2(hang(i), :), -Z2(i)); % 对每H6(i)行进行向右循环移位Z6（j）
    end 
end
g1 = bitxor(I2_All_1,D3_All_2);
countG = sum(g1(:) ~= 0);
%-------------------------------------------------------------------------%
%--------------------------------------对D3_All进行扩散逆过程--------------%
D3_All_1=zeros(size(D3_All_2));
D3_All_1(1)=mod(256*2+D3_All_2(1)-I2_All_256(1),256);
for i=2:3*M*N-A
    D3_All_1(i)=mod(256*2+D3_All_2(i)-I2_All_256(i)-D3_All_2(i-1),256);
end
h1 = bitxor(I2_All,D3_All_1);
countH = sum(h1(:) ~= 0);
%-------------------------------------------------------------------------%
%-----------------------------将矩阵D3_All_1重塑为D3_r,D3_g,D3_b-------------%
D3_r = D3_All_1(1:M*N-a);
D3_g = D3_All_1(M*N-a+1:2*M*N-a-b);
D3_b = D3_All_1(2*M*N-a-b+1:3*M*N-a-b-c);
D3_r = reshape(D3_r,[mr,nr]);
D3_g = reshape(D3_g,[mg,ng]);
D3_b = reshape(D3_b,[mb,nb]);
%-------------------------------------------------------------------------%
%------------------------对D3_r,D3_g,D3_b进行S_box_B解密得到D2_r,D2_g,D2_b-----------------%
S_box1=reshape(s_box_B_inversion(S_box_B,rowidx1,columnsidx1,rowidx11,columnsidx11),[16,16]);
S_box2=reshape(s_box_B_inversion(S_box_B,rowidx2,columnsidx2,rowidx22,columnsidx22),[16,16]);
S_box3=reshape(s_box_B_inversion(S_box_B,rowidx3,columnsidx3,rowidx33,columnsidx33),[16,16]);
S_box4=reshape(s_box_B_inversion(S_box_B,rowidx4,columnsidx4,rowidx44,columnsidx44),[16,16]);
S_box5=reshape(s_box_B_inversion(S_box_B,rowidx5,columnsidx5,rowidx55,columnsidx55),[16,16]);
S_box6=reshape(s_box_B_inversion(S_box_B,rowidx6,columnsidx6,rowidx66,columnsidx66),[16,16]);

D2_r = zeros(size(I1_rL));
D2_g = zeros(size(I1_gL));
D2_b = zeros(size(I1_bL));

    D3_rL=floor(D3_r/16)+1;
    D3_rH=mod(D3_r,16)+1;
    D3_gL=floor(D3_g/16)+1;
    D3_gH=mod(D3_g,16)+1;
    D3_bL=floor(D3_b/16)+1;
    D3_bH=mod(D3_b,16)+1;

for i=1:M*N-a
   switch(pr5(i))
       case 1
           D2_r(i) = S_box1(rowidx11(D3_rH(i)),columnsidx11(D3_rL(i)));
       case 2
           D2_r(i) = S_box2(rowidx22(D3_rH(i)),columnsidx22(D3_rL(i)));
       case 3
           D2_r(i) = S_box3(rowidx33(D3_rH(i)),columnsidx33(D3_rL(i)));
       case 4
           D2_r(i) = S_box4(rowidx44(D3_rH(i)),columnsidx44(D3_rL(i)));
       case 5
           D2_r(i) = S_box5(rowidx55(D3_rH(i)),columnsidx55(D3_rL(i)));
       case 6
           D2_r(i) = S_box6(rowidx66(D3_rH(i)),columnsidx66(D3_rL(i)));
   end
end
for i=1:M*N-b
   switch(pg5(i))
       case 1
           D2_g(i) = S_box1(rowidx11(D3_gH(i)),columnsidx11(D3_gL(i)));
       case 2
           D2_g(i) = S_box2(rowidx22(D3_gH(i)),columnsidx22(D3_gL(i)));
       case 3
           D2_g(i) = S_box3(rowidx33(D3_gH(i)),columnsidx33(D3_gL(i)));
       case 4
           D2_g(i) = S_box4(rowidx44(D3_gH(i)),columnsidx44(D3_gL(i)));
       case 5
           D2_g(i) = S_box5(rowidx55(D3_gH(i)),columnsidx55(D3_gL(i)));
       case 6
           D2_g(i) = S_box6(rowidx66(D3_gH(i)),columnsidx66(D3_gL(i)));
   end
end
for i=1:M*N-c
   switch(pb5(i))
       case 1
           D2_b(i) = S_box1(rowidx11(D3_bH(i)),columnsidx11(D3_bL(i)));
       case 2
           D2_b(i) = S_box2(rowidx22(D3_bH(i)),columnsidx22(D3_bL(i)));
       case 3
           D2_b(i) = S_box3(rowidx33(D3_bH(i)),columnsidx33(D3_bL(i)));
       case 4
           D2_b(i) = S_box4(rowidx44(D3_bH(i)),columnsidx44(D3_bL(i)));
       case 5
           D2_b(i) = S_box5(rowidx55(D3_bH(i)),columnsidx55(D3_bL(i)));
       case 6
           D2_b(i) = S_box6(rowidx66(D3_bH(i)),columnsidx66(D3_bL(i)));
   end
end

b1 = bitxor(I1_r,D2_r);
countB = sum(b1(:) ~= 0);
%-------------------------------------------------------------------------%
%---------------------------对S_sort_box进行S_box_A的解密------------------%
S_box_A_in = s_box_A_inversion(S_box_A);%先构造S_box_A的拟盒
    S_sort_box_L = floor(S_sort_box/16)+1;
    S_sort_box_H = mod(S_sort_box,16)+1;
    S_sort_back = zeros(1,A);
    for i=1:A          
              S_sort_back(i)=S_box_A_in(S_sort_box_H(i),S_sort_box_L(i));
    end
%-------------------------------------------------------------------------%
%---------------------------------对S_sort_back进行逆排序------------------%
[~, idx1] = sort(idx);% 
S_back = S_sort_back(idx1);%S_back=S
%-------------------------------------------------------------------------%
%---将a,b,c个元素分别放入D2_r,D2_g,D2_b中，得到[M,N]的矩阵D1_r,D1_g,D1_b----%
D1_r = zeros(M,N);
D1_g = zeros(M,N);
D1_b = zeros(M,N);

    [~,D1a_idx]=sort(Da);
    [~,D1b_idx]=sort(Db);
    [~,D1c_idx]=sort(Dc);
    
    S_back_r = S_back(1:a);
    S_back_g = S_back(a+1:a+b);
    S_back_b = S_back(a+b+1:a+b+c);
    
    e1 = bitxor(S_back_b,S_b);
    countE = sum(e1(:) ~= 0);
%     
    S_back_rr = S_back_r(D1a_idx);
    S_back_gg = S_back_g(D1b_idx);
    S_back_bb = S_back_b(D1c_idx);
    
for i=1:a+1 %将S_sort_back_sort中的a个元素插入到I3_r(Da(i))
    if(i==1)
        D1_r(1:D1a(i)-1)=D2_r(1:D1a(i)-1);
        D1_r(D1a(i))=S_back_rr(i);
    elseif(i<=a)
        D1_r(D1a(i-1)+1:D1a(i)-1)=D2_r(D1a(i-1)-i+2:D1a(i)-i);
        D1_r(D1a(i))=S_back_rr(i);
    else
        D1_r(D1a(a)+1:M*N)=D2_r(D1a(a)-a+1:mr*nr);
    end
end
for i=1:b+1 %将S_idx中的b个元素插入到I3_g(Db(i))
    if(i==1)
        D1_g(1:D1b(i)-1)=D2_g(1:D1b(i)-1);
        D1_g(D1b(i))=S_back_gg(i);
    elseif(i<=b)
        D1_g(D1b(i-1)+1:D1b(i)-1)=D2_g(D1b(i-1)-i+2:D1b(i)-i);
        D1_g(D1b(i))=S_back_gg(i);
    else
        D1_g(D1b(b)+1:M*N)=D2_g(D1b(b)-b+1:mg*ng);
    end
end
for i=1:c+1 %将S_idx中的c个元素插入到I3_b(Da(i))
    if(i==1)
        D1_b(1:D1c(i)-1)=D2_b(1:D1c(i)-1);
        D1_b(D1c(i))=S_back_bb(i);
    elseif(i<=c)
        D1_b(D1c(i-1)+1:D1c(i)-1)=D2_b(D1c(i-1)-i+2:D1c(i)-i);
        D1_b(D1c(i))=S_back_bb(i);
    else
        D1_b(D1c(c)+1:M*N)=D2_b(D1c(c)-c+1:mb*nb);
    end
end
D1(:,:,1)=uint8(D1_r);
D1(:,:,2)=uint8(D1_g);
D1(:,:,3)=uint8(D1_b);
figure(2)
imshow(D1);




