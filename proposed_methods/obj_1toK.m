function [sumK] = obj_1toK(B_est,ran,delta_f,valid_freq,c0)
% -----------------------------------------------------------------------------------------    
% The function is used to construct the objective funciton in the correlation-based scheme.
% version 1.0 - 01/07/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -----------------------------------------------------------------------------------------
sumK = 0;
for k = 1:size(B_est,2)
    ID = find(valid_freq(:,k)).';
    sumK = sumK + abs(normal_column(B_est(ID,k))'*normal_column(exp(-1j * 2 * pi * delta_f * (ID.'-1) * (2*ran/c0))));
end
end