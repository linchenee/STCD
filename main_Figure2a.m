clc
clear all
close all
addpath(genpath('./')) 

Mr = 10;         % number of receive antennas
N = 6;           % number of subcarriers
K = 6;           % number of OFDM symbols
p1 = 0;          % correspond to the parameter 'p_1' in the paper
p2 = 0.1:0.1:1;  % correspond to the parameter 'p_2' in the paper
num_trial = 500; % Monte Carlo trial number 
N1 = N*p1;
K1 = K*p1;

L1 = zeros(length(p2),num_trial); % generic upper bound of L in Theorem 1
L2 = zeros(length(p2),num_trial); % generic upper bound of L in Theorem 2
for k = 1:length(p2)
    Rate = p2(k);
    fprintf('%d\n',k);
    for trial = 1:num_trial
        rng(trial);
        mask = generate_mask_Type1(Rate,trial,Mr,N,K,N1,K1); % binary mask for the time-frequency resource
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
xlabel('|\Omega|/({\itN\itK})', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('Generic upper bound of \itL', 'FontName', 'Times New Roman','FontSize', 18);
legend1 = legend('Theorem 1','Theorem 2', 'FontName', 'Times New Roman','FontSize', 14);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6',...
    '0.7','0.8','0.9','1'},'XTick',[1,2,3,4,5,6,7,8,9,10],'XLim',[1,10],'YLim',[1,10],'YTick',[1,2,4,6,8,10]);
set(legend1,'Position',[0.150 0.796 0.237 0.114]);