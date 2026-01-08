function [InputMUSIC1, InputMUSIC2, InputMUSIC3] = prepare_input_for_MUSIC(Y_mask, mask, N1, K1)
% ------------------------------------------------------------------    
% The function is used to prepare the input for the MUSIC algorithm.
% version 1.0 - 01/06/2026
% Written by Lin Chen (lchen53@stevens.edu)
% ------------------------------------------------------------------
[Mr, N, K] = size(mask);
num_sample = sum(sum(mask(1,:,:)));

InputMUSIC1 = zeros(Mr,num_sample);
sample = 1;
for n = 1:N
    for k = 1:K 
        if mask(1,n,k) == 1
           InputMUSIC1(:,sample) = Y_mask(:,n,k);
           sample = sample + 1;
        end
    end
end

InputMUSIC2 = zeros(N1,Mr*K1);
sample = 1;
for m = 1:Mr
    for k = 1:K1 
       temp = Y_mask(m,1:N1,k).';
       InputMUSIC2(:,sample) = temp;
       sample = sample + 1;
    end
end

InputMUSIC3 = zeros(K1,Mr*N1);
sample = 1;
for m = 1:Mr
    for n = 1:N1 
       temp = squeeze(Y_mask(m,n,1:K1));
       InputMUSIC3(:,sample) = temp;
       sample = sample + 1;
    end
end

end