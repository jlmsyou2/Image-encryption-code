function inv_s_box = s_box_B_inversion (s_box,rowidx1,columnsidx1,rowidx11,columnsidx11)
%S_BOX_INVERSION  Invert S-box.
%
%   [INV_S_BOX] = S_BOX_INVERSION (S_BOX) 
%   creates the inverse S-box
%   from the previously created S-box.

%   Copyright 2001-2005, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

%   Version 1.0     30.05.2001


% Loop over all byte values
inv_s_box = zeros(16,16);
L=floor(s_box/16)+1;
H=mod(s_box,16)+1;
k=1;
for j = 1 : 16
    for i = 1:16
%         inv_s_box(H(i),L(i)) = rowidx1(mod(i,16))-1+(columnsidx1(floor(i/16)+1)-1)*16;
        inv_s_box(rowidx11(H(i+(j-1)*16)),columnsidx11(L(i+(j-1)*16))) = rowidx1(i)-1+(columnsidx1(j)-1)*16;
        k=k+1;
    end

    % Create the inverse S-box by taking the values 
    % of the elements of the S-Box as indices:
    % e.g.: s_box(00hex) = 63hex   ==>   inv_s_box(63hex) = 00hex
    % (except the fact, that Matlab vectors start at 1...)
%     inv_s_box(s_box(i) + 1) = i - 1;
    
end
a=1;
    


