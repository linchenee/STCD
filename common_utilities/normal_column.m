function [A,lambda] = normal_column(A)
    % ----------------------------------------------------------    
    % The function is used to normalize each column of a matrix.
    % version 1.0 - 01/07/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ----------------------------------------------------------
    lambda = sqrt(sum(abs(A).^2,1));
    M = size(A,1);
    A = A./(ones(M,1)*lambda);
end