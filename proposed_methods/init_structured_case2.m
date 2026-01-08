function [A_init,B_init,C_init] = init_structured_case2(Y_mask,mask,W,Input1,L,interval,minRV,maxR,maxV,gridNum,delta_f,lambda,Ts,c0)
% ---------------------------------------------------------------------------------------------------------------    
% The function is used for structured initialization of the STCD-2 method in Case-2, as described in Section V-C.
% version 1.0 - 01/07/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ---------------------------------------------------------------------------------------------------------------
[Mr,N,K] = size(Y_mask);
[Us,~,~] = svds(Input1,L);
[~,Z] = eig(pinv(Us(1:Mr-1,:))*Us(2:Mr,:));
vecZ = diag(Z).';  vecZ = vecZ./abs(vecZ);
Angle = asin(angle(vecZ)/pi)/pi*180;

gridR = linspace(minRV,maxR,gridNum); % range grids
gridV = linspace(minRV,maxV,gridNum); % velocity grids
A_init = exp(1j * pi * (0 : (Mr - 1)).' * sin(Angle * pi / 180));
pinv_A_init = pinv(A_init);

est_ran = zeros(1,L); % estimated range parameters
est_vel = zeros(1,L); % estimated velocity parameters
Ball = exp(-1j * 2 * pi * delta_f * (0:N-1).' * (2 * gridR / c0));
Dall = exp(1j * 2 * pi * (0:K-1).' * (2 * gridV / lambda) *Ts);
%% correlation-based scheme
for i = 1:L
    corr = 0;
    B_est = zeros(N,K);
    valid_freq = zeros(N,K);
    for k=1:K
        valid_freq(:,k) = double(mask(1,:,k)==1).';
        ID = find(valid_freq(:,k)).';
        B_est(:,k) = (pinv_A_init(i,:)*squeeze(Y_mask(:,:,k))).';
        corr = corr + abs(normal_column(B_est(ID,k))'*normal_column(Ball(ID,:)));
    end
    [~, index] = max(corr);
    est_ran(1,i) = gridR(index); % coarse estimate

    corr = 0;
    C_est = zeros(K,N);
    valid_time = zeros(K,N);
    for n=1:N
        valid_time(:,n) = double(squeeze(mask(1,n,:))==1);
        ID = find(valid_time(:,n)).';
        C_est(:,n) = (pinv_A_init(i,:)*squeeze(Y_mask(:,n,:))).';
        Ctemp = diag(W(:,ID).'*A_init(:,i))*Dall(ID,:);
        corr = corr + abs(normal_column(C_est(ID,n))'*normal_column(Ctemp));
    end
    [~, index] = max(corr);
    est_vel(1,i) = gridV(index); % coarse estimate
    
    %% refine the range estimate and velocity estimate
    obj_func = @(ran) 1/obj_1toK(B_est,ran,delta_f,valid_freq,c0);
    est_ran(1,i) = fminsearchbnd(obj_func,est_ran(i),max(est_ran(i)-interval,minRV),min(est_ran(i)+interval,maxR));
    obj_func = @(vel) 1/obj_1toN(C_est,vel,lambda,valid_time,Ts,A_init(:,i),W);
    est_vel(1,i) = fminsearchbnd(obj_func,est_vel(i),max(est_vel(i)-interval,minRV),min(est_vel(i)+interval,maxV));
end
B_init = exp(-1j * 2 * pi * delta_f * (0 : N - 1).' * (2 * est_ran / c0));
C_init = zeros(K,L);
OMEGA = find(mask);
Z1 = zeros(length(OMEGA),L);
for i = 1 : L
    D = diag(exp(1j*2*pi*(0:K-1)*(2 * est_vel(i) / lambda)*Ts));
    C_init(:, i) = D * W.'*A_init(:,i);
    temp = permute(reshape(kr(A_init(:,i),B_init(:,i))*C_init(:,i).',N,Mr,K),[2,1,3]);
    Z1(:,i) = temp(OMEGA);
end
Gain_init = pinv(Z1)*Y_mask(OMEGA); 
C_init = C_init*diag(Gain_init);
end
