function [Angle,Range,Veloc,Doppl,Delay,Gain] = target_parameter(num_tar,k,lambda,c0)
    % ---------------------------------------------------    
    % The function is used to generate target parameters.
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ---------------------------------------------------
    rng(k);
    Angle = 90 * rand(num_tar, 1) - 45;          % It corresponds to the parameter 'theta_l' in the paper.
    Range = 200 * rand(num_tar, 1);              % It corresponds to the parameter 'R_l' in the paper.
    Veloc = 30 * rand(num_tar, 1);               % It corresponds to the parameter 'v_l' in the paper.
    Doppl = 2 * Veloc / lambda;                  % It corresponds to the parameter 'nu_l' in the paper.
    Delay = 2 * Range.' / c0;                    % It corresponds to the parameter 'tau_l' in the paper.
    power = rand(num_tar,1).*3;
    phase = rand(num_tar,1)*360;
    Gain = 10.^(power/10).*exp(1i*phase*pi/180); % It corresponds to the parameter 'alpha_l' in the paper.
end