function [mask] = generate_mask_Type2(p2,k,Mr,N,K,N1,K1)
    % ------------------------------------------------------------------------    
    % The function is used to generate the Type-2 mask, as shown in Fig. 3(b).
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ------------------------------------------------------------------------
    while 1
        rng(k);
        %96/12=8;
        %56/14=4;
        mask_temp = (rand(8,4) < p2);
        mask = zeros(Mr,N, K);
        for i = 1:8
            for j = 1:4
                if mask_temp(i,j) == 1
                   mask(:,(i-1)*12+1:i*12,(j-1)*14+1:j*14) = 1;
                end
            end
        end
        mask(:,1:N1,1:K1) = 1;
        if min(sum(squeeze(mask(1,:,:)),2)) > 0
           break;
        else 
           k = k + 300;
        end
    end
    % imagesc(squeeze(mask(1,:,:)));
    % figure();
    % imagesc(squeeze(mask(2,:,:)));
end