clc
clear all
close all
addpath(genpath('./')) 

global c0 N K lambda delta_f Ts Mr Mt

system_parameter_loc

%% Scenario setting
SNRs = [0,5,10,15,20,25,30];
p1 = 0; % correspond to the parameter 'p_1' in the paper
p2 = 0.6; % correspond to the parameter 'p_2' in the paper
N1 = N*p1;
K1 = K*p1;
num_trial = 1000; % Monte Carlo trial number 

%% Parameter setting of the STCD/UTCD method 
I = 20;     % iteration number in the alternating minimization (AM)
J = 2;      % iteration number in the projected gradient descent (PGD)
xi = 1e-8;  % termination tolerance in the AM
eta = 1e-7; % stepsize in the PGD
thr = 1;    % termination tolerance in the PGD

E = zeros(numel(SNRs),num_trial); 
RMSE_range0 = E; RMSE_range2 = E; RMSE_range3 = E; RMSE_range4 = E; RMSE_range2part2 = E; rcrb_range = E;
RMSE_angle0 = E; RMSE_angle2 = E; RMSE_angle3 = E; RMSE_angle4 = E; RMSE_angle2part2 = E; rcrb_angle = E;
RMSE_veloc0 = E; RMSE_veloc2 = E; RMSE_veloc3 = E; RMSE_veloc4 = E; RMSE_veloc2part2 = E; rcrb_veloc = E;

minRV = 0;  % minimum of range and velocity
maxR = 200; % maximum of range
maxV = 30;  % maximum of velocity
gridNum = 1000; % number of girds 
interval = 0.5;

W = dftmtx(K); W = W(1:Mt, :); % precoder
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

        %% MUSIC method (Irregular)
        [Input1,~,~] = prepare_input_for_MUSIC(Y_mask, mask, N1, K1);
        Angle4 = MUSICstep1(Input1, size(Input1,2), L); % angle estimate
        RMSE_angle4(i_snr, trial) = cal_RMSE(Angle4*pi/180, Angle*pi/180);

        %% STCD-2 method (Irregular)
        [A_init,B_init,C_init] = init_structured_case2(Y_mask,mask,W,Input1,L,interval,minRV,maxR,maxV,gridNum,delta_f,lambda,Ts,c0);
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2,Angle2,Veloc2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle2(i_snr, trial) = cal_RMSE(Angle2*pi/180, Angle*pi/180);
        RMSE_range2(i_snr, trial) = cal_RMSE(Range2, Range);
        RMSE_veloc2(i_snr, trial) = cal_RMSE(Veloc2, Veloc);

        %% STCD-1 method (Irregular)
        [A_init,B_init,C_init] = init_unstructured(Y_mask,L,find(mask),xi);
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2part2,Angle2part2,Veloc2part2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle2part2(i_snr, trial) = cal_RMSE(Angle2part2*pi/180, Angle*pi/180);
        RMSE_range2part2(i_snr, trial) = cal_RMSE(Range2part2, Range);
        RMSE_veloc2part2(i_snr, trial) = cal_RMSE(Veloc2part2, Veloc);

        %% UTCD method (Irregular)
        [A3, B3, C3] = UTCD(Y_mask,find(mask),I,xi,A_init,B_init,C_init);
        [Range3,Angle3,Veloc3] = para_extra_for_TCD(A3,B3,C3,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle3(i_snr, trial) = cal_RMSE(Angle3*pi/180, Angle*pi/180);
        RMSE_range3(i_snr, trial) = cal_RMSE(Range3, Range);
        RMSE_veloc3(i_snr, trial) = cal_RMSE(Veloc3, Veloc); 

        [rcrb_angle(i_snr,trial),rcrb_range(i_snr,trial),rcrb_veloc(i_snr,trial)] = ...
            cal_CRB(Angle.',Delay,Doppl.',Gain.',Mr,N,K,sigma_2,mask,delta_f,Ts,W,L,c0,lambda);
    end
end

close all
thr = 0.05;
[RMSE_ang0,~] = mean_RMSE_success_rate(RMSE_angle0,thr,num_trial);
[RMSE_ang2,~] = mean_RMSE_success_rate(RMSE_angle2,thr,num_trial);
[RMSE_ang3,~] = mean_RMSE_success_rate(RMSE_angle3,thr,num_trial);
[RMSE_ang4,~] = mean_RMSE_success_rate(RMSE_angle4,thr,num_trial);
[RMSE_ang2part2,~] = mean_RMSE_success_rate(RMSE_angle2part2,thr,num_trial);

thr = 1;
[RMSE_ran0,~] = mean_RMSE_success_rate(RMSE_range0,thr,num_trial);
[RMSE_ran2,~] = mean_RMSE_success_rate(RMSE_range2,thr,num_trial);
[RMSE_ran3,~] = mean_RMSE_success_rate(RMSE_range3,thr,num_trial);
[RMSE_ran2part2,~] = mean_RMSE_success_rate(RMSE_range2part2,thr,num_trial);

thr = 2;
[RMSE_vel0,~] = mean_RMSE_success_rate(RMSE_veloc0,thr,num_trial);
[RMSE_vel2,~] = mean_RMSE_success_rate(RMSE_veloc2,thr,num_trial);
[RMSE_vel3,~] = mean_RMSE_success_rate(RMSE_veloc3,thr,num_trial);
[RMSE_vel2part2,~] = mean_RMSE_success_rate(RMSE_veloc2part2,thr,num_trial);

figure1 = figure(1);
semilogy(squeeze(mean(RMSE_ang4,2)),'-v', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on 
semilogy(squeeze(mean(RMSE_ang3,2)),'k-o', 'LineWidth', 3); hold on 
semilogy(squeeze(mean(RMSE_ang2part2,2)),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on 
semilogy(squeeze(mean(RMSE_ang2,2)),'r-s', 'LineWidth', 3); hold on 
semilogy(squeeze(mean(rcrb_angle,2)),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (rad)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0','5','10','15','20','25','30'},'YLim',[0.00002,0.015]);
legend('MUSIC (Irregular)', 'UTCD (Irregular)', 'STCD-1 (Irregular)', 'STCD-2 (Irregular)', 'RCRB (Irregular)', 'FontName', 'Times New Roman', 'FontSize', 14);
set(legend,'Position',[0.677249998552459 0.717261901582991 0.212500002895083 0.28095238731021],'FontSize',15);

figure2 = figure(2);
semilogy(squeeze(mean(RMSE_ran3,2)),'k-o', 'LineWidth', 3); hold on 
semilogy(squeeze(mean(RMSE_ran2part2,2)),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on
semilogy(squeeze(mean(RMSE_ran2,2)),'r-s', 'LineWidth', 3); hold on 
semilogy(squeeze(mean(rcrb_range,2)),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on 
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0','5','10','15','20','25','30'}, 'YLim',[0.000657112953401,0.1]); 
legend('UTCD (Irregular)', 'STCD-1 (Irregular)', 'STCD-2 (Irregular)', 'RCRB (Irregular)', 'FontName', 'Times New Roman', 'FontSize', 14);
set(legend,'Position',[0.684 0.766 0.205 0.219],'FontSize',15);

figure3 = figure(3);
semilogy(squeeze(mean(RMSE_vel3,2)),'k-o', 'LineWidth', 3); hold on
semilogy(squeeze(mean(RMSE_vel2part2,2)),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on 
semilogy(squeeze(mean(RMSE_vel2,2)),'r-s', 'LineWidth', 3); hold on
semilogy(squeeze(mean(rcrb_veloc,2)),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on 
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m/s)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0','5','10','15','20','25','30'}, 'YLim',[0.005184522989455,0.42]);
legend('UTCD (Irregular)', 'STCD-1 (Irregular)', 'STCD-2 (Irregular)', 'RCRB (Irregular)', 'FontName', 'Times New Roman', 'FontSize', 14);
set(legend,'Position',[0.677249998552459 0.770523806980678 0.212500002895083 0.226428576514835],'FontSize',15);

