function [sr] = root_MUSIC(Y,k)
% ------------------
% Root-MUSIC method
% ------------------
[L, samples] = size(Y);
[eigvec, ~] = eig(Y*Y'/samples);

%% Computing coefficients for polynomial construction
En = eigvec(:,1:L-k); % noise eigenspace
A = En*En';
cf = zeros(L-1,1);
for i=1:L-1
    cf(i) = sum(diag(A,i));
end
% Complete set of coefficients
CF = [flipud(cf); sum(diag(A)); conj(cf)];

% Computing roots 
rts = roots(CF);
rts(find(abs(rts)>1)) = 0;
[~,ind] = maxk(abs(rts),k);
true_rts = rts(ind);
%% From roots to transmitted code
sr = angle(true_rts);
end