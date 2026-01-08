function [Y,sigma_2] = add_noise_to_echo(X,trial,SNR,N,K,Mr)
    % -----------------------------------------------------------------    
    % The function is used to add the noise Z to the echo X in Eq. (5).
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % -----------------------------------------------------------------
    rng(trial);
    Noise = permute((randn(N,K,Mr)+1j*randn(N,K,Mr))/sqrt(2),[3,1,2]);
    sigma = (norm(X(:),'fro')/norm(Noise(:),'fro'))*10^(-SNR/20);
    sigma_2 = sigma^2;
    Z = sigma*Noise;
    Y = X + Z;
end