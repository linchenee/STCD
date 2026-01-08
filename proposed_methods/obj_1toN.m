function [sumN] = obj_1toN(C_est,vel,lambda,valid_time,Ts,A_initi,W)
% -----------------------------------------------------------------------------------------    
% The function is used to construct the objective funciton in the correlation-based scheme.
% version 1.0 - 01/07/2026
% Written by Lin Chen (lchen53@stevens.edu)
% -----------------------------------------------------------------------------------------
sumN = 0;
for n = 1:size(C_est,2)
    ID = find(valid_time(:,n)).';
    sumN = sumN + abs(normal_column(C_est(ID,n))'*normal_column( (W(:,ID).'*A_initi).*exp(1j*2*pi*(ID.'-1)* (2*vel/ lambda) *Ts)    ));
end
end