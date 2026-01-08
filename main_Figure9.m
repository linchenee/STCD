clc
clear all
close all
addpath(genpath('./')) 

global c0 N K lambda delta_f Ts Mr Mt

system_parameter_loc

%% Scenario setting
SNR = 5;
p1 = 0.5; % correspond to the parameter 'p_1' in the paper
p2 = [0.1 0.3 0.5 0.7 0.9]; % correspond to the parameter 'p_2' in the paper
N1 = N*p1;
K1 = K*p1;
num_trial = 1000; % Monte Carlo trial number 

%% Parameter setting of the STCD/UTCD method 
I = 20;     % iteration number in the alternating minimization (AM)
J = 2;      % iteration number in the projected gradient descent (PGD)
xi = 1e-8;  % termination tolerance in the AM
eta = 1e-7; % stepsize in the PGD
thr = 2;    % termination tolerance in the PGD

E = zeros(numel(p2),num_trial); 
RMSE_range0 = E; RMSE_range1 = E; RMSE_range2 = E; RMSE_range3 = E; RMSE_range4 = E; RMSE_range2part2 = E; rcrb_range = E;
RMSE_angle0 = E; RMSE_angle1 = E; RMSE_angle2 = E; RMSE_angle3 = E; RMSE_angle4 = E; RMSE_angle2part2 = E; rcrb_angle = E;
RMSE_veloc0 = E; RMSE_veloc1 = E; RMSE_veloc2 = E; RMSE_veloc3 = E; RMSE_veloc4 = E; RMSE_veloc2part2 = E; rcrb_veloc = E;

W = dftmtx(K); W = W(1:Mt, :); % precoder
for i_p2 = 1:length(p2)
    Rate = p2(i_p2);
    fprintf('The %dth p2\n',i_p2);
    parfor trial = 1:num_trial
        % fprintf('The %dth trial\n',trial);
        rng(trial);
        L = randi(4);
        [Angle,Range,Veloc,Doppl,Delay,Gain] = target_parameter(L,trial,lambda,c0); % ground-truth target parameters 
        X = generate_echo(Angle,Delay,Doppl,Gain,W,L,N,K,delta_f,Ts,Mr,Mt); % echo tensor in Eq. (6)
        mask = generate_mask_Type1(Rate,trial,Mr,N,K,N1,K1); % binary mask for the time-frequency resource
        [Y, sigma_2] = add_noise_to_echo(X,trial,SNR,N,K,Mr); % Y = X+Z in Eq. (5)
        Y_mask = Y .* mask; % observation

        %% MUSIC method
        [Input1,Input2,~] = prepare_input_for_MUSIC(Y_mask, mask, N1, K1);
        Angle4 = MUSICstep1(Input1, size(Input1,2), L); % angle estimate
        Range4 = MUSICstep2(Input2, size(Input2,2), L, delta_f, c0); % range estimate
        RMSE_angle4(i_p2, trial) = cal_RMSE(Angle4*pi/180, Angle*pi/180);
        RMSE_range4(i_p2, trial) = cal_RMSE(Range4, Range);

        %% VCPD method (Regular)
        [delay1, angle1] = VCPD(permute(Y(:,1:N1,1:K1),[2,1,3]),L);
        [Range1,Angle1,Veloc1] = para_extra_for_VCPD_IMDF(angle1, delay1, Y(:,1:N1,1:K1),L,W(:,1:K1),c0,lambda,delta_f,Ts,Mt,Mr,N1,K1);
        Doppl1 = 2 * Veloc1 / lambda; % Doppler estimate
        RMSE_angle1(i_p2, trial) = cal_RMSE(Angle1*pi/180, Angle*pi/180);
        RMSE_range1(i_p2, trial) = cal_RMSE(Range1, Range);
        RMSE_veloc1(i_p2, trial) = cal_RMSE(Veloc1, Veloc);

        %% IMDF method (Regular)
        [wx, wy] = IMDF_multi3D(Y(:,1:N1,1:K1), L);
        [Range0,Angle0,Veloc0] = para_extra_for_VCPD_IMDF(wx, wy, Y(:,1:N1,1:K1),L,W(:,1:K1),c0,lambda,delta_f,Ts,Mt,Mr,N1,K1);
        RMSE_angle0(i_p2, trial) = cal_RMSE(Angle0*pi/180, Angle*pi/180);
        RMSE_range0(i_p2, trial) = cal_RMSE(Range0, Range);
        RMSE_veloc0(i_p2, trial) = cal_RMSE(Veloc0, Veloc);

        %% STCD-2 method (Irregular)
        [A_init,B_init,C_init] = init_structured_case1(Y,W,L,Mr,N,K,N1,K1,Angle1,Range1,Doppl1,delta_f,Ts,c0);
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2,Angle2,Veloc2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle2(i_p2, trial) = cal_RMSE(Angle2*pi/180, Angle*pi/180);
        RMSE_range2(i_p2, trial) = cal_RMSE(Range2, Range);
        RMSE_veloc2(i_p2, trial) = cal_RMSE(Veloc2, Veloc);

        %% STCD-1 method (Irregular)
        [A_init,B_init,C_init] = init_unstructured(Y_mask,L,find(mask),xi);
        [A2, B2, C2] = STCD(Y_mask,L,find(mask),I,J,xi,eta,thr,A_init,B_init,C_init);
        [Range2part2,Angle2part2,Veloc2part2] = para_extra_for_TCD(A2,B2,C2,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle2part2(i_p2, trial) = cal_RMSE(Angle2part2*pi/180, Angle*pi/180);
        RMSE_range2part2(i_p2, trial) = cal_RMSE(Range2part2, Range);
        RMSE_veloc2part2(i_p2, trial) = cal_RMSE(Veloc2part2, Veloc);

        %% UTCD method (Irregular)
        [A3, B3, C3] = UTCD(Y_mask,find(mask),I,xi,A_init,B_init,C_init);
        [Range3,Angle3,Veloc3] = para_extra_for_TCD(A3,B3,C3,W,L,c0,lambda,K,delta_f,Ts,Mt);
        RMSE_angle3(i_p2, trial) = cal_RMSE(Angle3*pi/180, Angle*pi/180);
        RMSE_range3(i_p2, trial) = cal_RMSE(Range3, Range);
        RMSE_veloc3(i_p2, trial) = cal_RMSE(Veloc3, Veloc);

        [rcrb_angle(i_p2,trial),rcrb_range(i_p2,trial),rcrb_veloc(i_p2,trial)] = ...
            cal_CRB(Angle.',Delay,Doppl.',Gain.',Mr,N,K,sigma_2,mask,delta_f,Ts,W,L,c0,lambda);
    end
end

close all
thr = 0.05;
[RMSE_ang0,~] = mean_RMSE_success_rate(RMSE_angle0,thr,num_trial);
[RMSE_ang1,~] = mean_RMSE_success_rate(RMSE_angle1,thr,num_trial);
[RMSE_ang2,~] = mean_RMSE_success_rate(RMSE_angle2,thr,num_trial);
[RMSE_ang3,~] = mean_RMSE_success_rate(RMSE_angle3,thr,num_trial);
[RMSE_ang4,~] = mean_RMSE_success_rate(RMSE_angle4,thr,num_trial);
[RMSE_ang2part2,~] = mean_RMSE_success_rate(RMSE_angle2part2,thr,num_trial);

thr = 1;
[RMSE_ran0,~] = mean_RMSE_success_rate(RMSE_range0,thr,num_trial);
[RMSE_ran1,~] = mean_RMSE_success_rate(RMSE_range1,thr,num_trial);
[RMSE_ran2,~] = mean_RMSE_success_rate(RMSE_range2,thr,num_trial);
[RMSE_ran3,~] = mean_RMSE_success_rate(RMSE_range3,thr,num_trial);
[RMSE_ran4,~] = mean_RMSE_success_rate(RMSE_range4,thr,num_trial);
[RMSE_ran2part2,~] = mean_RMSE_success_rate(RMSE_range2part2,thr,num_trial);

thr = 2;
[RMSE_vel0,~] = mean_RMSE_success_rate(RMSE_veloc0,thr,num_trial);
[RMSE_vel1,~] = mean_RMSE_success_rate(RMSE_veloc1,thr,num_trial);
[RMSE_vel2,~] = mean_RMSE_success_rate(RMSE_veloc2,thr,num_trial);
[RMSE_vel3,~] = mean_RMSE_success_rate(RMSE_veloc3,thr,num_trial);
[RMSE_vel2part2,~] = mean_RMSE_success_rate(RMSE_veloc2part2,thr,num_trial);

figure1=figure(1);
x1=semilogy(squeeze(mean(RMSE_ang4,2)),'-v', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on 
x2=semilogy(squeeze(mean(RMSE_ang1,2)),'b:>', 'LineWidth', 3); hold on 
x3=semilogy(squeeze(mean(RMSE_ang0,2)),':^', 'LineWidth', 3, 'Color', [0.72,0.27,1.00]); hold on 
x4=semilogy(squeeze(mean(RMSE_ang3,2)),'k-o', 'LineWidth', 3); hold on
x5=semilogy(squeeze(mean(RMSE_ang2part2,2)),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on 
x6=semilogy(squeeze(mean(RMSE_ang2,2)),'r-s', 'LineWidth', 3); hold on 
x7=semilogy(squeeze(mean(rcrb_angle,2)),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('{\itp}_2', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (rad)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0.1','0.3','0.5','0.7','0.9'},'XLim',[1,5],'YLim',[0.00032,0.01],'YTick',[1e-3,1e-2]); 
empty_patch = patch([0 0 1 1], [0 1 1 0], 'w', 'EdgeColor', 'none'); hold on;
legend([x1, x2, x3, x4], 'MUSIC (Irregular)', 'VCPD (Regular)', 'IMDF (Regular)', 'UTCD (Irregular)', 'FontName', 'Times New Roman', 'FontSize', 14, 'EdgeColor', 'white', ...
    'Position',[0.297261904761906 0.776428566523963 0.343928577457155 0.218809528714135]);
temp = axes('position', get(gca, 'position'), 'visible', 'off');
legend(temp, [x5, x6, x7, empty_patch], 'STCD-1 (Irregular)', 'STCD-2 (Irregular)', 'RCRB (Irregular)', ' ', 'FontName', 'Times New Roman', 'FontSize', 14, 'EdgeColor', 'white', ...
    'Position',[0.641883810092229 0.778174598269995 0.351071434770311 0.218809528714135]);
annotation(figure1,'textbox',[0.296785714285714 0.774761904761905 0.695 0.220952380952384],'FitBoxToText','off');

figure2 = figure(2);
x1=semilogy(squeeze(mean(RMSE_ran4,2)),':v', 'LineWidth', 3, 'Color', [0.93,0.69,0.13]); hold on 
x2=semilogy(squeeze(mean(RMSE_ran1,2)),'b:>', 'LineWidth', 3); hold on 
x3=semilogy(squeeze(mean(RMSE_ran0,2)),':^', 'LineWidth', 3, 'Color', [0.72,0.27,1.00]); hold on 
x4=semilogy(squeeze(mean(RMSE_ran3,2)),'k-o', 'LineWidth', 3); hold on
x5=semilogy(squeeze(mean(RMSE_ran2part2,2)),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on 
x6=semilogy(squeeze(mean(RMSE_ran2,2)),'r-s', 'LineWidth', 3); hold on
x7=semilogy(squeeze(mean(rcrb_range,2)),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on 
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('{\itp}_2', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0.2','0.3','0.5','0.7','0.9'},'XLim',[1,5],'YLim',[0.0094,0.33]); 
empty_patch = patch([0 0 1 1], [0 1 1 0], 'w', 'EdgeColor', 'none'); hold on;
legend([x1, x2, x3, x4], 'MUSIC (Regular)', 'VCPD (Regular)', 'IMDF (Regular)', 'UTCD (Irregular)', 'FontName', 'Times New Roman', 'FontSize', 14, 'EdgeColor', 'white', ...
    'Position',[0.297261904761906 0.776428566523963 0.343928577457155 0.218809528714135]);
temp = axes('position', get(gca, 'position'), 'visible', 'off');
legend(temp, [x5, x6, x7, empty_patch], 'STCD-1 (Irregular)', 'STCD-2 (Irregular)', 'RCRB (Irregular)', ' ', 'FontName', 'Times New Roman', 'FontSize', 14, 'EdgeColor', 'white', ...
    'Position',[0.641883810092229 0.778174598269995 0.351071434770311 0.218809528714135]);
annotation(figure2,'textbox',[0.296785714285714 0.774761904761905 0.695 0.220952380952384],'FitBoxToText','off');

figure3 = figure(3);
x1=semilogy(squeeze(mean(RMSE_vel1,2)),'b:>', 'LineWidth', 3); hold on
x2=semilogy(squeeze(mean(RMSE_vel0,2)),':^', 'LineWidth', 3, 'Color', [0.72,0.27,1.00]); hold on 
x3=semilogy(squeeze(mean(RMSE_vel3,2)),'k-o', 'LineWidth', 3); hold on 
x4=semilogy(squeeze(mean(RMSE_vel2part2,2)),'r-+', 'LineWidth', 3, 'Color', [1 0.412 0.161]); hold on 
x5=semilogy(squeeze(mean(RMSE_vel2,2)),'r-s', 'LineWidth', 3); hold on
x6=semilogy(squeeze(mean(rcrb_veloc,2)),'-d', 'LineWidth', 3, 'Color', [0.39,0.83,0.07]); hold on 
grid on
set(gca, 'GridColor', [0, 0, 0], 'GridLineWidth', 1);
xlabel('{\itp}_2', 'FontName', 'Times New Roman','FontSize', 18);
ylabel('RMSE (m/s)', 'FontName', 'Times New Roman','FontSize', 18);
set(gca, 'FontName', 'Times New Roman','FontSize', 18, 'XTickLabel',{'0.1','0.3','0.5','0.7','0.9'},'XLim',[1,5],'YLim',[0.072,1]);
legend([x1, x2, x3],'VCPD (Regular)', 'IMDF (Regular)', 'UTCD (Irregular)', 'FontName', 'Times New Roman', 'FontSize', 14, 'EdgeColor', 'white',...
    'Position',[0.310238089550109 0.82793650425805 0.329642862830843 0.166190479868934]);
temp = axes('position', get(gca, 'position'), 'visible', 'off');
legend(temp, [x4, x5, x6], 'STCD-1 (Irregular)', 'STCD-2 (Irregular)', 'RCRB (Irregular)', 'FontName', 'Times New Roman', 'FontSize', 14, 'EdgeColor', 'white',...
    'Position',[0.639523803324927 0.828888885210431 0.351071434770312 0.166190479868934]);
annotation(figure3,'textbox',[0.311357142857143 0.826190476190476 0.679 0.16761904761905],'FitBoxToText','off');