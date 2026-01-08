function [para1_esti, para2_esti] = VCPD(X,R)
% ---------------------------------------------------------------------------------------------    
% Vandermonde CANDECOMP/PARAFAC decomposition (VCPD) method in [3].
% [3] R. Zhang, L. Cheng, S. Wang, Y. Lou, Y. Gao, W. Wu, and D. W. K. Ng, "Integrated sensing 
%     and communication with massive MIMO: A unified tensor approach for channel and target 
%     parameter estimation," IEEE Trans. Wireless Commun., vol. 23, no. 8, pp. 8571–8587, 2024.
% version 1.0 - 01/06/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ---------------------------------------------------------------------------------------------
I1 = size(X,1);          
I2 = size(X,2);
I3 = size(X,3);
[K1,L1,K2,L2] = findPairs_2d(I1, I2, I3);
X3 = reshape(permute(X,[2,1,3]),I1*I2,I3);
Smooth = spatial_smoothing(X3, L1, L2, K1, K2);

[U,~,~] = svds(Smooth,R); 
U1 = U(1:(K1-1)*K2,:);
U2 = U(K2+1:K1*K2,:);
[M1,Z1] = eig(pinv(U1)*U2);
vecZ1 = (diag(Z1)).';
vecZ1 = vecZ1./abs(vecZ1);
para1_esti = angle(vecZ1);
A1_esti = zeros(I1,R);
for i = 1:I1
  A1_esti(i,:) = vecZ1.^(i-1);
end

vecZ2 = zeros(1,R);
AK2_esti = zeros(K2,R);
for i = 1:R
  AK2_esti(:,i) = kron((A1_esti(1:K1,i))',eye(K2))*U*M1(:,i);
  vecZ2(1,i) = pinv(AK2_esti(1:K2-1,i))*AK2_esti(2:K2,i);
end
vecZ2 = vecZ2./abs(vecZ2);
para2_esti = angle(vecZ2);
end

function J = spatial_smoothing(X, L1, L2, K1, K2)
  J=[];
  for i=1:L1
      for j=1:L2
         J=[J, kron([zeros(K1,i-1),eye(K1),zeros(K1,(L1-i))], [zeros(K2,j-1),eye(K2),zeros(K2,(L2-j))])*X];
      end
  end
end

function [K1,L1,K2,L2] = findPairs_2d(I1,I2,I3)
%% find the optimal K1, L1, K2, and L2
    r=0;
    for k1=1:I1
         l1=I1-k1+1;
         for k2=1:I2
             l2=I2-k2+1;
             temp = min( (k1-1)*k2,l1*l2*I3 );
             if temp>r
               r=temp;
               K1=k1; L1=l1; K2=k2; L2=l2;
             end
         end
    end
end