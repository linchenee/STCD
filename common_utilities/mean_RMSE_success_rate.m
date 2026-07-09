function [mean_RMSE,success_rate] = mean_RMSE_success_rate(RMSEs,thr,num_trial)
    temp = sort(RMSEs,2);
    success_rate = mean(temp<thr,2);
    mean_RMSE = sum(temp.*(temp<thr),2)./(success_rate*num_trial);
end