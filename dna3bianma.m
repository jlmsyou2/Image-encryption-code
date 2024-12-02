function [s] = dna3bianma(a,c)
[M,N]=size(a);
a=reshape(a,M*N,1);
a=dec2bin(a,8);
a=a-48;
% a=reshape(a,8,M*N);
% s=zeros(4,M*N);
s=zeros(M,N*4);
n=1;
for i=1:M*N
    if c(i)==1
        for j=1:2:8
           if a(i,j)==0&&a(i,j+1)==1
                s(n)='A';
            end
            if a(i,j)==1&&a(i,j+1)==0
                s(n)='T';
            end
            if a(i,j)==0&&a(i,j+1)==0
                s(n)='G';
            end
            if a(i,j)==1&&a(i,j+1)==1
                s(n)='C';
            end 
            n=n+1;
        end
    end
    if c(i)==2
        for j=1:2:8
           if a(i,j)==1&&a(i,j+1)==1
                s(n)='A';
            end
            if a(i,j)==0&&a(i,j+1)==0
                s(n)='T';
            end
           if a(i,j)==1&&a(i,j+1)==0
                s(n)='G';
           end
           if a(i,j)==0&&a(i,j+1)==1
                s(n)='C';
           end
           n=n+1;
        end
    end
    if c(i)==3
        for j=1:2:8
           if a(i,j)==1&&a(i,j+1)==1
                s(n)='A';
            end
            if a(i,j)==0&&a(i,j+1)==0
                s(n)='T';
            end
           if a(i,j)==0&&a(i,j+1)==1
                s(n)='G';
           end
           if a(i,j)==1&&a(i,j+1)==0
                s(n)='C';
           end
           n=n+1;
        end
    end
    if c(i)==4
        for j=1:2:8
           if a(i,j)==1&&a(i,j+1)==0
                s(n)='A';
            end
            if a(i,j)==0&&a(i,j+1)==1
                s(n)='T';
            end
           if a(i,j)==0&&a(i,j+1)==0
                s(n)='G';
           end
           if a(i,j)==1&&a(i,j+1)==1
                s(n)='C';
           end
           n=n+1;
        end
    end
    if c(i)==5
        for j=1:2:8
           if a(i,j)==1&&a(i,j+1)==0
                s(n)='A';
            end
            if a(i,j)==0&&a(i,j+1)==1
                s(n)='T';
            end
           if a(i,j)==1&&a(i,j+1)==1
                s(n)='G';
           end
           if a(i,j)==0&&a(i,j+1)==0
                s(n)='C';
           end
           n=n+1;
        end
    end
    if c(i)==6
        for j=1:2:8
           if a(i,j)==0&&a(i,j+1)==0
                s(n)='A';
            end
            if a(i,j)==1&&a(i,j+1)==1
                s(n)='T';
            end
           if a(i,j)==0&&a(i,j+1)==1
                s(n)='G';
           end
           if a(i,j)==1&&a(i,j+1)==0
                s(n)='C';
           end
           n=n+1;
        end
    end
    if c(i)==7
        for j=1:2:8
           if a(i,j)==0&&a(i,j+1)==1
                s(n)='A';
            end
            if a(i,j)==1&&a(i,j+1)==0
                s(n)='T';
            end
           if a(i,j)==1&&a(i,j+1)==1
                s(n)='G';
           end
           if a(i,j)==0&&a(i,j+1)==0
                s(n)='C';
           end
           n=n+1;
        end
    end
    if c(i)==8
        for j=1:2:8
           if a(i,j)==0&&a(i,j+1)==0
                s(n)='A';
            end
            if a(i,j)==1&&a(i,j+1)==1
                s(n)='T';
            end
           if a(i,j)==1&&a(i,j+1)==0
                s(n)='G';
           end
           if a(i,j)==0&&a(i,j+1)==1
                s(n)='C';
           end
           n=n+1;
        end
    end
end
s=char(s);
end
