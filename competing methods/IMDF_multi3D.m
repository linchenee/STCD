%% If you find this code useful, please cite our paper:
% C. Qian, X. Fu, N. D. Sidiropoulos and Y. Yang, "Tensor-Based Channel 
% Estimation for Dual-Polarized Massive MIMO Systems," in IEEE Transactions 
% on Signal Processing. doi: 10.1109/TSP.2018.2873506

%% Thanks.

function [wx, wy] = IMDF_multi3D(X, P)

[N1, N2, N3] = size(X);
[K1,~,K2,~] = findPairs_2d(N1, N2);

Ztot = [];
for i = 1:N3
    X1 = squeeze(X(:,:,i));
    X1s = IMDFS2(X1, N1, N2, K1, K2);
    [LL1,LL2]=size(X1s);
    x1s = X1s(:);
    Y1s = reshape(conj(x1s(end:-1:1)), LL1, LL2);
    Z = [X1s, Y1s];
    Ztot = [Ztot, Z];
end

[U, ~, ~] = svd(Ztot);
us = U(:, 1:P);
J1 = kron(eye(K1-1,K1), eye(K2));
J2 = kron([zeros(K1-1,1),eye(K1-1)], eye(K2));

U1 = J1*us;
U2 = J2*us;
[Tsp, ~] = eig((U1'*U1)\U1'*U2);
Ahat = us*Tsp;

A11 = J1*Ahat.'';
A12 = J2*Ahat.'';

for i = 1:P
    wx(i) = -angle(A11(:,i)'*A12(:,i));
end

J1 = kron(eye(K1), eye(K2-1,K2));
J2 = kron(eye(K1), [zeros(K2-1,1),eye(K2-1)]);
A21 = J1*Ahat.'';
A22 = J2*Ahat.'';
for i = 1:P
    wy(i) = -angle(A21(:,i)'*A22(:,i));
end

% Ax = exp(1j*[0:Mt.left-1]'*wx(:)');
% Ay = exp(1j*[0:Mt.right-1]'*wy(:)');
% At = krb(Ay.'',Ax.'');
    
end

function Z=IMDFS2(X, N1, N2, K1, K2)
L1=N1-K1+1;
L2=N2-K2+1;
for i1=1:L1
    X1(:,:,i1)=X(i1:N1-L1+i1,:);
end
for i2=1:L2
    Y(:,:,:,i2)=X1(:,i2:N2-L2+i2,:);
end
Y2=permute(Y, [2,1,3,4]);
Z=reshape(Y2, K1*K2, L1*L2);
end

function [k1,k2,l1,l2] = findPairs_2d(M,N)
%% find the optimal k1 and l1

Fold = 0;
for k1 = floor(M/2):M
    k2 = M + 1 - k1;
    for l1 = floor(N/2):N
        l2 = N + 1 - l1;
        F = 0;
        if M == 1
            while min(2*l2,l1-1)>=F || min(l1*2,l2-1)>=F
                F = F + 1;
            end
        elseif M == 2
            while (2*k2*l2>=F && (k1-1)*l1>=F )
                F = F + 1;
            end
        else
            while (8*k2*l2>=F && k1*(l1-1)>=F )
                F = F + 1;
            end
        end
        if F>Fold
            l1new = l1;
            k1new = k1;
            Fold = F - 1;
        end
    end
end
l1 = l1new;
k1 = k1new;
k2 = M + 1 - k1;
l2 = N + 1 - l1;
end
