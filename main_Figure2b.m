clc
clear all
close all
addpath(genpath('./')) 

Mr = 10;         % number of receive antennas
K_vec = 2:18;    % number of OFDM symbols
p1 = 0;          % correspond to the parameter 'p_1' in the paper
p2 = 0.3;        % correspond to the parameter 'p_2' in the paper
num_trial = 500; % Monte Carlo trial number
N1 = 0;          % (N1=N*p1)
K1 = 0;          % (K1=K*p1)

L1 = zeros(length(K_vec),num_trial); % generic upper bound of L in Theorem 1
L2 = zeros(length(K_vec),num_trial); % generic upper bound of L in Theorem 2
for k = 1:length(K_vec)
    K = K_vec(k);
    N = K; % number of subcarriers
    fprintf('%d\n',k);
    for trial = 1:num_trial
        rng(trial);
        mask = generate_mask_Type1(p2,trial,Mr,N,K,N1,K1);
        L1(k,trial) = CP_rank_bound1(mask);
        for i = 1:Mr-1
            for j = 1:Mr-1
                if Mr+2-i-j >= 2
                   L2(k,trial) = max(L2(k,trial), CP_rank_bound2(mask,i,j));
                end
            end
        end
    end
end

figure(1);
plot(mean(L1,2),'MarkerSize',8,'LineWidth',2,'Marker','square'); hold on
plot(mean(L2,2),'MarkerSize',8,'LineWidth',2,'Marker','o'); hold on
xlabel('\itK', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('Generic upper bound of \itL', 'FontName', 'Times New Roman','FontSize', 18);
legend1 = legend('Theorem 1','Theorem 2', 'FontName', 'Times New Roman','FontSize', 14);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'2','4','6','8','10','12','14','16','18'},...
    'XTick',[1,3,5,7,9,11,13,15,17],'XLim',[1,17],'YLim',[1,10],'YTick',[1,2,4,6,8,10]);
set(legend1,'Position',[0.150 0.797 0.237 0.114]);
