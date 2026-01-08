function y = cal_RMSE(vec_a, vec_b)
    % ------------------------------------------
    % The function is used to calculate RMSE.
    % version 1.0 - 01/06/2026
    % Written by Lin Chen (lchen53@stevens.edu)
    % ------------------------------------------
    vec_a = vec_a(:);
    vec_b = vec_b(:);
    mat_a = perms(vec_a).';
    mat_b = repmat(vec_b, 1, size(mat_a, 2));
    y = min(sqrt((1/numel(vec_a)) * sum(abs(mat_a - mat_b).^2, 1)));
end