function [Y] = PGD(X,W11,W1,p,r,step_size,option,temp1,temp2)
% For comprehensive documentation of this function, visit:
% https://github.com/linchenee/ES-CPD/blob/main/functions/utilities/PGD.m

[n1,n2] = size(X);

Z = Hankel_transform( X-2*step_size*(W11.*(X*temp1-temp2)) ); % Equation (18)

% Perform the rank-r truncated tensor SVD in Equation (19) and 
% the inverse Hankel transform in Euqation (20) jointly
Y = zeros(n1,n2);
for k = 1:n2
    [u,s,v] = svds(Z(:,:,k),r,'largest',option);
    temp = ifft(fft(u,p,1).*fft(conj(v),p,1),[],1);
    Y(:,k) = sum(W1.*temp(1:n1,:)*s,2);
end

end