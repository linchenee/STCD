clc
clear all
close all
addpath(genpath('./'))

global c0 N K lambda delta_f Ts Mr Mt

system_parameter_loc

%% Scenario setting
SNRs = [-40,-30,-20,-10,0,10,20,30];
p1 = 0.5; % correspond to the parameter 'p_1' in the paper
p2 = 0.5; % correspond to the parameter 'p_2' in the paper
L = 1; % single target
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
RMSE_range0 = E; RMSE_range1 = E; RMSE_range2 = E; RMSE_range3 = E; RMSE_range4 = E; rcrb_range = E; 
RMSE_range2part2 = E; RMSE_range2part3 = E; RMSE_range2part4 = E; RMSE_range3part2 = E; rcrb_rangepart2 = E;

RMSE_angle0 = E; RMSE_angle1 = E; RMSE_angle2 = E; RMSE_angle3 = E; RMSE_angle4 = E; rcrb_angle = E;
RMSE_angle2part2 = E; RMSE_angle3part2 = E;  RMSE_angle2part3 = E; RMSE_angle2part4 = E; rcrb_anglepart2 = E;

RMSE_veloc0 = E; RMSE_veloc1 = E; RMSE_veloc2 = E; RMSE_veloc3 = E; RMSE_veloc4 = E; rcrb_veloc = E;
RMSE_veloc2part2 = E; RMSE_veloc2part3 = E; RMSE_veloc2part4 = E; RMSE_veloc3part2 = E; rcrb_velocpart2 = E;

for i_snr = 1:length(SNRs)
    SNR = SNRs(i_snr);
    fprintf('The %dth SNR\n',i_snr);
    parfor trial = 1:num_trial
        % fprintf('The %dth trial\n',trial);
        rng(trial);
        W = exp(1j*randn(Mt,K)); % precoder          
        [Angle,Range,Veloc,Doppl,Delay,Gain] = target_parameter(L,trial,lambda,c0); % ground-truth target parameters   
        X = generate_echo(Angle,Delay,Doppl,Gain,W,L,N,K,delta_f,Ts,Mr,Mt); % echo tensor in Eq. (6)
        mask = generate_mask_Type1(p2,trial,Mr,N,K,N1,K1); % binary mask for the time-frequency resource
        [Y, sigma_2] = add_noise_to_echo(X,trial,SNR,N,K,Mr); % Y = X+Z in Eq. (5)
        Y_mask = Y .* mask; % observation

        %% MUSIC method
        [Input1,Input2,~] = prepare_input_for_MUSIC(Y_mask, mask, N1, K1);
        Angle4 = asind((1/(2*pi))*2*root_MUSIC(Input1,L)); % angle estimate
        Range4 = para_extra_for_MUSIC(root_MUSIC(Input2,L),delta_f,c0); % range estimate
        RMSE_angle4(i_snr, trial) = cal_RMSE(Angle4*pi/180, Angle*pi/180); % calculate RMSE
        RMSE_range4(i_snr, trial) = cal_RMSE(Range4, Range);

        %% VCPD method (Regular)
        [delay1, angle1] = VCPD(permute(Y(:,1:N1,1:K1),[2,1,3]),L);
        [Range1,Angle1,Veloc1] = para_extra_for_VCPD_IMDF(angle1, delay1, Y(:,1:N1,1:K1),L,W(:,1:K1),c0,lambda,delta_f,Ts,Mt,Mr,N1,K1);
        Doppl1 = 2 * Veloc1 / lambda; % Doppler estimate
        RMSE_angle1(i_snr, trial) = cal_RMSE(Angle1*pi/180, Angle*pi/180);
        RMSE_range1(i_snr, trial) = cal_RMSE(Range1, Range);
        RMSE_veloc1(i_snr, trial) = cal_RMSE(Veloc1, Veloc);

        %% IMDF method (Regular)
        [wx, wy] = IMDF_multi3D(Y(:,1:N1,1:K1), L);
        [Range0,Angle0,Veloc0] = para_extra_for_VCPD_IMDF(wx, wy, Y(:,1:N1,1:K1),L,W(:,1:K1),c0,lambda,delta_f,Ts,Mt,Mr,N1,K1);
        RMSE_angle0(i_snr, trial) = cal_RMSE(Angle0*pi/180, Angle*pi/180);
        RMSE_range0(i_snr, trial) = cal_RMSE(Range0, Range);
        RMSE_veloc0(i_snr, trial) = cal_RMSE(Veloc0, Veloc);

        %% STCD-2 method (Irregular)
        [A_init,B_init,C_init] = init_structured_case1(Y,W,L,Mr,N,K,N1,K1,Angle1,Range1,Doppl1,delta_f,Ts,c0);
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2,Angle2,Veloc2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle2(i_snr, trial) = cal_RMSE(Angle2*pi/180, Angle*pi/180);
        RMSE_range2(i_snr, trial) = cal_RMSE(Range2, Range);
        RMSE_veloc2(i_snr, trial) = cal_RMSE(Veloc2, Veloc);

        %% STCD-2 method (Regular)
        B_init = B_init(1:N1,:);
        C_init = C_init(1:K1,:);
        [A2, B2, C2] = STCD(Y_mask(:,1:N1,1:K1),L,find(mask(:,1:N1,1:K1)),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2part2,Angle2part2,Veloc2part2] = para_extra_for_TCD(A2,B2,C2,W(:,1:K1),L,c0,lambda,K1,delta_f,Ts,Mt);
        RMSE_angle2part2(i_snr, trial) = cal_RMSE(Angle2part2*pi/180, Angle*pi/180);
        RMSE_range2part2(i_snr, trial) = cal_RMSE(Range2part2, Range);
        RMSE_veloc2part2(i_snr, trial) = cal_RMSE(Veloc2part2, Veloc);

        %% STCD-1 method (Irregular)
        [A_init,B_init,C_init] = init_unstructured(Y_mask,L,find(mask),xi);
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2part3,Angle2part3,Veloc2part3] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle2part3(i_snr, trial) = cal_RMSE(Angle2part3*pi/180, Angle*pi/180);
        RMSE_range2part3(i_snr, trial) = cal_RMSE(Range2part3, Range);
        RMSE_veloc2part3(i_snr, trial) = cal_RMSE(Veloc2part3, Veloc);

        %% UTCD method (Irregular)
        [A3, B3, C3] = UTCD(Y_mask,find(mask),I,xi,A_init,B_init,C_init);
        [Range3,Angle3,Veloc3] = para_extra_for_TCD(A3,B3,C3,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle3(i_snr, trial) = cal_RMSE(Angle3*pi/180, Angle*pi/180);
        RMSE_range3(i_snr, trial) = cal_RMSE(Range3, Range);
        RMSE_veloc3(i_snr, trial) = cal_RMSE(Veloc3, Veloc);

        %% STCD-1 method (Regular)
        [A_init,B_init,C_init] = init_unstructured(Y_mask(:,1:N1,1:K1),L,find(mask(:,1:N1,1:K1)),xi);
        [A2, B2, C2] = STCD(Y_mask(:,1:N1,1:K1),L,find(mask(:,1:N1,1:K1)),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2part4,Angle2part4,Veloc2part4] = para_extra_for_TCD(A2,B2,C2,W(:,1:K1),L,c0,lambda,K1,delta_f,Ts,Mt);
        RMSE_angle2part4(i_snr, trial) = cal_RMSE(Angle2part4*pi/180, Angle*pi/180);
        RMSE_range2part4(i_snr, trial) = cal_RMSE(Range2part4, Range);
        RMSE_veloc2part4(i_snr, trial) = cal_RMSE(Veloc2part4, Veloc);

        %% UTCD method (Regular)
        [A3, B3, C3] = UTCD(Y_mask(:,1:N1,1:K1),find(mask(:,1:N1,1:K1)),I,xi,A_init,B_init,C_init);
        [Range3part2,Angle3part2,Veloc3part2] = para_extra_for_TCD(A3,B3,C3,W(:,1:K1),L,c0,lambda,K1,delta_f,Ts,Mt);
        RMSE_angle3part2(i_snr, trial) = cal_RMSE(Angle3part2*pi/180, Angle*pi/180);
        RMSE_range3part2(i_snr, trial) = cal_RMSE(Range3part2, Range);
        RMSE_veloc3part2(i_snr, trial) = cal_RMSE(Veloc3part2, Veloc);

        %% RCRB (Irregular)
        [rcrb_angle(i_snr,trial),rcrb_range(i_snr,trial),rcrb_veloc(i_snr,trial)] = ...
            cal_CRB(Angle.',Delay,Doppl.',Gain.',Mr,N,K,sigma_2,mask,delta_f,Ts,W,L,c0,lambda);
 
        %% RCRB (Regular)
        [rcrb_anglepart2(i_snr,trial),rcrb_rangepart2(i_snr,trial),rcrb_velocpart2(i_snr,trial)] = ...
            cal_CRB(Angle.',Delay,Doppl.',Gain.',Mr,N1,K1,sigma_2,mask(:,1:N1,1:K1),delta_f,Ts,W(:,1:K1),L,c0,lambda);
    end
end

close all
x = (-40:10:30)';
figure(1);
semilogy(x,mean(RMSE_angle4,2),'-v', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on 
semilogy(x,mean(RMSE_angle1,2),'b:>', 'LineWidth', 3); hold on 
semilogy(x,mean(RMSE_angle0,2),':^', 'LineWidth', 3, 'Color', [0.72,0.27,1.00]); hold on 
semilogy(x,mean(RMSE_angle3part2,2),'k:o', 'LineWidth', 3); hold on 
semilogy(x,mean(RMSE_angle3,2),'k-o', 'LineWidth', 3); hold on 
semilogy(x,mean(RMSE_angle2part4,2),'r:+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on
semilogy(x,mean(RMSE_angle2part3,2),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on
semilogy(x,mean(RMSE_angle2part2,2),'r:s', 'LineWidth', 3); hold on 
semilogy(x,mean(RMSE_angle2,2),'r-s', 'LineWidth', 3); hold on 
semilogy(x,mean(rcrb_anglepart2,2),':d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on 
semilogy(x,mean(rcrb_angle,2),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on 
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (rad)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'-40','-30','-20','-10','0','10','20','30'},'XTick',[-40,-30,-20,-10,0,10,20,30],'XLim',[-40,30],'YLim',[1.05e-5,1.3]);
legend1 = legend('MUSIC (Irregular)','VCPD (Regular)','IMDF (Regular)','UTCD (Regular)','UTCD (Irregular)','STCD-1 (Regular)','STCD-1 (Irregular)','STCD-2 (Regular)','STCD-2 (Irregular)','RCRB (Regular)','RCRB (Irregular)', 'FontName', 'Times New Roman','FontSize', 14);
set(legend1,'Position',[0.6092 0.4103 0.35107 0.58714]);
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);

figure(2);
semilogy(x,mean(RMSE_range4,2),':v', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on
semilogy(x,mean(RMSE_range1,2),'b:>', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_range0,2),':^', 'LineWidth', 3, 'Color', [0.72,0.27,1.00]); hold on
semilogy(x,mean(RMSE_range3part2,2),'k:o', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_range3,2),'k-o', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_range2part4,2),'r:+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on
semilogy(x,mean(RMSE_range2part3,2),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on
semilogy(x,mean(RMSE_range2part2,2),'r:s', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_range2,2),'r-s', 'LineWidth', 3); hold on
semilogy(x,mean(rcrb_rangepart2,2),':d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on
semilogy(x,mean(rcrb_range,2),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'-40','-30','-20','-10','0','10','20','30'},'XTick',[-40,-30,-20,-10,0,10,20,30],'XLim',[-40,30],'YLim',[3e-4,5e2]);
legend1 = legend('MUSIC (Regular)','VCPD (Regular)','IMDF (Regular)','UTCD (Regular)','UTCD (Irregular)','STCD-1 (Regular)','STCD-1 (Irregular)','STCD-2 (Regular)','STCD-2 (Irregular)','RCRB (Regular)','RCRB (Irregular)', 'FontName', 'Times New Roman','FontSize', 14);
set(legend1,'Position',[0.6070 0.4103 0.35107 0.58714]);
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);

figure(3);
semilogy(x,mean(RMSE_veloc1,2),'b:>', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_veloc0,2),':^', 'LineWidth', 3, 'Color', [0.72,0.27,1.00]); hold on
semilogy(x,mean(RMSE_veloc3part2,2),'k:o', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_veloc3,2),'k-o', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_veloc2part4,2),'r:+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on
semilogy(x,mean(RMSE_veloc2part3,2),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on
semilogy(x,mean(RMSE_veloc2part2,2),'r:s', 'LineWidth', 3); hold on
semilogy(x,mean(RMSE_veloc2,2),'r-s', 'LineWidth', 3); hold on
semilogy(x,mean(rcrb_velocpart2,2),':d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on
semilogy(x,mean(rcrb_veloc,2),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on
xlabel('SNR (dB)', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m/s)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'-40','-30','-20','-10','0','10','20','30'},'XTick',[-40,-30,-20,-10,0,10,20,30],'XLim',[-40,30],'YLim',[1.52e-3,8e2]);
legend1 = legend('VCPD (Regular)','IMDF (Regular)','UTCD (Regular)','UTCD (Irregular)','STCD-1 (Regular)','STCD-1 (Irregular)','STCD-2 (Regular)','STCD-2 (Irregular)','RCRB (Regular)','RCRB (Irregular)', 'FontName', 'Times New Roman','FontSize', 14);
set(legend1,'Position',[0.6 0.4610 0.3511 0.5345]);
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
