clc
clear all
close all
addpath(genpath('./'))

global c0 N K lambda delta_f Ts Mr Mt

system_parameter_loc

%% Scenario setting
SNRs = [0,5,10,15,20,25,30];
p1 = 0.5; % correspond to the parameter 'p_1' in the paper
p2 = 0.5; % correspond to the parameter 'p_2' in the paper
N1 = N*p1;
K1 = K*p1;
num_trial = 1000; % Monte Carlo trial number 

%% Parameter setting of the STCD/UTCD method 
I = 20;     % iteration number in the alternating minimization (AM)
xi = 1e-8;  % termination tolerance in the AM
eta = 1e-7; % stepsize in the PGD
thr = 2;    % termination tolerance in the PGD

E = zeros(numel(SNRs),num_trial); 
RMSE_range0 = E; RMSE_range1 = E; RMSE_range2 = E; RMSE_range3 = E;
RMSE_angle0 = E; RMSE_angle1 = E; RMSE_angle2 = E; RMSE_angle3 = E;
RMSE_veloc0 = E; RMSE_veloc1 = E; RMSE_veloc2 = E; RMSE_veloc3 = E;

W = dftmtx(K); W = W(1:Mt, :);
for i_snr = 1:length(SNRs)
    SNR = SNRs(i_snr);
    fprintf('The %dth SNR\n',i_snr);
    parfor trial = 1:num_trial
        % fprintf('The %dth trial\n',trial);
        rng(trial);
        L = randi(4);
        [Angle,Range,Veloc,Doppl,Delay,Gain] = target_parameter(L,trial,lambda,c0); % ground-truth target parameters 
        X = generate_echo(Angle,Delay,Doppl,Gain,W,L,N,K,delta_f,Ts,Mr,Mt); % echo tensor in Eq. (6)
        mask = generate_mask_Type1(p2,trial,Mr,N,K,N1,K1); % binary mask for the time-frequency resource
        [Y, sigma_2] = add_noise_to_echo(X,trial,SNR,N,K,Mr); % Y = X+Z in Eq. (5)
        Y_mask = Y .* mask; % observation

        %% VCPD method (Regular)
        [delay1, angle1] = VCPD(permute(Y(:,1:N1,1:K1),[2,1,3]),L);
        [Range1,Angle1,Veloc1] = para_extra_for_VCPD_IMDF(angle1, delay1, Y(:,1:N1,1:K1),L,W(:,1:K1),c0,lambda,delta_f,Ts,Mt,Mr,N1,K1);
        Doppl1 = 2 * Veloc1 / lambda; % Doppler estimate

        %% STCD-2 method (Irregular)
        [A_init,B_init,C_init] = init_structured_case1(Y,W,L,Mr,N,K,N1,K1,Angle1,Range1,Doppl1,delta_f,Ts,c0);

        J = 1;
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2,Angle2,Veloc2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle0(i_snr, trial) = cal_RMSE(Angle2*pi/180, Angle*pi/180);
        RMSE_range0(i_snr, trial) = cal_RMSE(Range2, Range);
        RMSE_veloc0(i_snr, trial) = cal_RMSE(Veloc2, Veloc);

        J = 2;
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2,Angle2,Veloc2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle1(i_snr, trial) = cal_RMSE(Angle2*pi/180, Angle*pi/180);
        RMSE_range1(i_snr, trial) = cal_RMSE(Range2, Range);
        RMSE_veloc1(i_snr, trial) = cal_RMSE(Veloc2, Veloc);

        J = 3;
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2,Angle2,Veloc2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle2(i_snr, trial) = cal_RMSE(Angle2*pi/180, Angle*pi/180);
        RMSE_range2(i_snr, trial) = cal_RMSE(Range2, Range);
        RMSE_veloc2(i_snr, trial) = cal_RMSE(Veloc2, Veloc);

        J = 4;
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2,Angle2,Veloc2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle3(i_snr, trial) = cal_RMSE(Angle2*pi/180, Angle*pi/180);
        RMSE_range3(i_snr, trial) = cal_RMSE(Range2, Range);
        RMSE_veloc3(i_snr, trial) = cal_RMSE(Veloc2, Veloc);
    end
end

close all
thr = 0.05;
[RMSE_ang0,~] = mean_RMSE_success_rate(RMSE_angle0,thr,num_trial);
[RMSE_ang1,~] = mean_RMSE_success_rate(RMSE_angle1,thr,num_trial);
[RMSE_ang2,~] = mean_RMSE_success_rate(RMSE_angle2,thr,num_trial);
[RMSE_ang3,~] = mean_RMSE_success_rate(RMSE_angle3,thr,num_trial);

thr = 1;
[RMSE_ran0,~] = mean_RMSE_success_rate(RMSE_range0,thr,num_trial);
[RMSE_ran1,~] = mean_RMSE_success_rate(RMSE_range1,thr,num_trial);
[RMSE_ran2,~] = mean_RMSE_success_rate(RMSE_range2,thr,num_trial);
[RMSE_ran3,~] = mean_RMSE_success_rate(RMSE_range3,thr,num_trial);

thr = 2;
[RMSE_vel0,~] = mean_RMSE_success_rate(RMSE_veloc0,thr,num_trial);
[RMSE_vel1,~] = mean_RMSE_success_rate(RMSE_veloc1,thr,num_trial);
[RMSE_vel2,~] = mean_RMSE_success_rate(RMSE_veloc2,thr,num_trial);
[RMSE_vel3,~] = mean_RMSE_success_rate(RMSE_veloc3,thr,num_trial);

figure(1);
semilogy(squeeze(mean(RMSE_ang0,2)),'-s', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on
semilogy(squeeze(mean(RMSE_ang1,2)),'b-^', 'LineWidth', 3); hold on 
semilogy(squeeze(mean(RMSE_ang2,2)),'r-+', 'LineWidth', 3); hold on 
semilogy(squeeze(mean(RMSE_ang3,2)),'k-o', 'LineWidth', 3); hold on 
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (rad)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0','5','10','15','20','25','30'},'YLim',[0.00002,0.0017],'YTick',[1e-4,1e-3]); 
legend('{\itJ}=1','{\itJ}=2','{\itJ}=3','{\itJ}=4', 'FontName', 'Times New Roman','FontSize', 14);

figure(2);
semilogy(squeeze(mean(RMSE_ran0,2)),'-s', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on 
semilogy(squeeze(mean(RMSE_ran1,2)),'b-^', 'LineWidth', 3); hold on 
semilogy(squeeze(mean(RMSE_ran2,2)),'r-+', 'LineWidth', 3); hold on
semilogy(squeeze(mean(RMSE_ran3,2)),'k-o', 'LineWidth', 3); hold on
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0','5','10','15','20','25','30'},'YLim',[0.001,0.07],'YTick',[1e-3,1e-2]); 
legend('{\itJ}=1','{\itJ}=2','{\itJ}=3','{\itJ}=4', 'FontName', 'Times New Roman','FontSize', 14);

figure(3);
semilogy(squeeze(mean(RMSE_vel0,2)),'-s', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on
semilogy(squeeze(mean(RMSE_vel1,2)),'b-^', 'LineWidth', 3); hold on
semilogy(squeeze(mean(RMSE_vel2,2)),'r-+', 'LineWidth', 3); hold on
semilogy(squeeze(mean(RMSE_vel3,2)),'k-o', 'LineWidth', 3); hold on
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m/s)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0','5','10','15','20','25','30'},'YLim',[0.006,0.4],'YTick',[1e-2,1e-1]);
legend('{\itJ}=1','{\itJ}=2','{\itJ}=3','{\itJ}=4', 'FontName', 'Times New Roman','FontSize', 14);