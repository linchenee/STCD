function [Range] = para_extra_for_MUSIC(input,delta_f,c0)
    % -------------------------------------------------------------------------------    
    % The function is used to extract target's range parameter for the MUSIC method.
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % -------------------------------------------------------------------------------
    Range = zeros(length(input),1); 
    for i = 1 : length(input)
        if input(i) > 0
            Range(i,1) = (c0/(2*delta_f)) * ((2*pi - input(i))/(2*pi));
        else
            Range(i,1) = (c0/(2*delta_f)) * (-input(i)/(2*pi));
        end
    end
end