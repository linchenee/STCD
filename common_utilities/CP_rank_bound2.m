function [result] = CP_rank_bound2(mask,V1,V2)
    % ------------------------------------------------------------------------
    % The function is used to calculate the generic CP rank bound in Theorem 2
    % version 1.0 - 01/08/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ------------------------------------------------------------------------
    [dim1,dim2,dim3] = size(mask); % dim1=Mr, dim2=N, dim3=K
    X_mode3 = reshape(permute(mask,[2,1,3]),dim1*dim2,dim3);

    U2 = dim1+2-V2-V1;
    U1 = U2 + V2 - 1;

    Smooth1 = spatial_smoothing(X_mode3, U1, V1, dim2).';
    Smooth2 = spatial_smoothing(Smooth1, U2, V2, dim3).';
    X_new = permute(reshape(Smooth2,V1*dim2,U2,V2*dim3),[2,1,3]);

    Dim2 = V1*dim2;
    Dim3 = V2*dim3;
    dsel = 0; 
    for j1=1:Dim2-1
        for j2=j1+1:Dim2
            for k1=1:Dim3-1
                for k2=k1+1:Dim3
                    if X_new(1,j1,k1)==1 && X_new(1,j1,k2) && X_new(1,j2,k1) && X_new(1,j2,k2)
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
    result = min(L,U2);
end    

function J = spatial_smoothing(X, L3, K3, n2)
    J = [];
    for i=1:L3
         J=[J, kron([zeros(K3,i-1),eye(K3),zeros(K3,(L3-i))], eye(n2))*X];
    end
end