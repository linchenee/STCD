function [T] = construct_T(n1,n2)
    % ------------------------------------------
    % Construct the matrix T in Eq. (20)
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ------------------------------------------
    if mod(n1,2)==0
     T_temp = repmat([1:n1/2,n1/2:-1:1]',1,n2);
    else
     T_temp = repmat([1:(n1+1)/2,(n1-1)/2:-1:1]',1,n2);
    end
    T = ones(n1,n2)./T_temp;
end