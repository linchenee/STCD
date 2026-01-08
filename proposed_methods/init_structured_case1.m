function [A_init,B_init,C_init] = init_structured_case1(Y,W,L,Mr,N,K,N1,K1,Angle,Range,Doppl,delta_f,Ts,c0)
    % ---------------------------------------------------------------------------------------------------------------    
    % The function is used for structured initialization of the STCD-2 method in Case-1, as described in Section V-C.
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ---------------------------------------------------------------------------------------------------------------
    A_init = exp(1j * pi * (0 : (Mr - 1)).' * sin(Angle * pi / 180));
    B_init = exp(-1j * 2 * pi * delta_f * (0 : N - 1).' * (2 * Range / c0));
    C_init = zeros(K,L);
    S1 = zeros(Mr*N1*K1,L);
    for i = 1 : L
        D = diag(exp(1j*2*pi*(0:K-1)*Doppl(i)*Ts));
        C_init(:, i) = D * W.'*A_init(:,i);
        S1(:,i) = vec(permute(reshape(kr(A_init(:,i),B_init(1:N1,i))*C_init(1:K1,i).',N1,Mr,K1),[2,1,3]));
    end
    Gain_init = pinv(S1)*vec(Y(:,1:N1,1:K1)); 
    C_init = C_init*diag(Gain_init);
end