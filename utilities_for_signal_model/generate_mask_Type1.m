function [mask] = generate_mask_Type1(p2,trial,Mr,N,K,N1,K1)
    % ------------------------------------------------------------------------    
    % The function is used to generate the Type-1 mask, as shown in Fig. 3(a).
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ------------------------------------------------------------------------
    rng(trial);
    mask_1 = (rand(N, K) < p2);
    mask = zeros(Mr, N, K);
    for anten = 1:Mr
     mask(anten,:,:) = mask_1;
     mask(anten,1:N1,1:K1) = 1;
    end
end    