function [result] = CP_rank_bound1(mask)
% ------------------------------------------------------------------------
% The function is used to calculate the generic CP rank bound in Theorem 1
% version 1.0 - 01/08/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------------------------------------
I1 = size(mask,1); % Mr     
I2 = size(mask,2); % N
I3 = size(mask,3); % K
dsel = 0; 
for j1=1:I2-1
    for j2=j1+1:I2
        for k1=1:I3-1
            for k2=k1+1:I3
                if mask(1,j1,k1)==1 && mask(1,j1,k2)==1 && mask(1,j2,k1)==1 && mask(1,j2,k2)==1
                   dsel = dsel+1;
                end
            end
        end
    end
end
L = 1;
while 1
    if dsel >= L*(L-1)/2
       L = L + 1;
    else
       break;
    end
end
L = L - 1;
result = min(I1,L);
end