%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to detect MOMs (Metastable Oscillatory Modes) in different 
% frequency bands for a selected point in the parameter space and calculate
% Shannon Entropy of the LFP phase covariance matrices throughout time 
%  
%
%  Needs: bandpasshopf.m, conver_back_to_time.m, subplot_tight.m
%
% July 2021
% Joana Cabral, Francesca Castaldo 
% joanacabral@med.uminho.pt
% francesca.castaldo.20@ucl.ac.uk
% 
% Adapted in July 2024
% Marta Xavier
% marta.xavier@tecnico.ulisboa.pt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 

% Frequency bands to analyze (in Hz)
delta = [0.5 4];
theta = [4 8];
alpha = [8 13];
beta = [13 30];
low_pass = 30;

folder = 'C:\Users\marta\OneDrive\Documentos\LASEEB\Simulations\FEP_MOMs\';

K='1E1'; MD='10';
simulation_file = {strcat('Hopf_Simu_K', K, '_MD0.mat') strcat('Hopf_Simu_K', K, '_MD', MD, '.mat')};

folder_pars = fullfile(folder, strcat('K', K, '_MD', MD)); mkdir(folder_pars);
folder_img = fullfile(folder_pars, 'Imgs'); mkdir(folder_img);
    

%% 1 - Read and define baseline 
% Load data in 0 delay case, intermediate coupling to define the baseline

% NOTES - BASELINE
% If MD (mean delay) is 0, units keep oscillating at their natural frequency
% Since some areas are more coupled together than others, even with 0 MD
% these areas may exhibit more power across frequencies that is purely due
% to noisy interactions - hence a different threshold is defined for each node and band

% Zsave describes the state of each of the N oscillators at each moment in 
% time (with resolution dt_save); the real part corresponds to the LFP

load(fullfile(folder, simulation_file{1})); % a = -5; MD = 0

N = size(Zsave, 1);
Order = [1:2:N N:-2:2];
Zsave  = Zsave./(5*std(Zsave(:)));
Zsave = Zsave(Order, :);
Zdelta = zeros(size(Zsave));
Ztheta = zeros(size(Zsave));
Zalpha = zeros(size(Zsave));
Zbeta  = zeros(size(Zsave));
Zfilt  = zeros(size(Zsave));

% Bandpass filter (via FFT) to obtain power at each frequency band 
% Cycle through nodes 
for n = 1 : N
    
    % Syntax is bandpasshopf( data, [f_min f_max], fs ) 
    Zdelta(n, :) = bandpasshopf(Zsave(n, :), delta , 1/dt_save);
    Ztheta(n, :) = bandpasshopf(Zsave(n, :), theta , 1/dt_save);
    Zalpha(n, :) = bandpasshopf(Zsave(n, :), alpha , 1/dt_save);
    Zbeta(n, :)  = bandpasshopf(Zsave(n, :), beta  , 1/dt_save);
    Zfilt(n, :)=   bandpasshopf(Zsave(n, :),[0.01 low_pass], 1/dt_save);
    
end % nodes 
clear Zsave

% Remove the first and last second after bandpass filtering 
Zdelta(:, [1:1/dt_save end-1/dt_save:end]) = [];
Ztheta(:, [1:1/dt_save end-1/dt_save:end]) = [];
Zalpha(:, [1:1/dt_save end-1/dt_save:end]) = [];
Zbeta(:, [1:1/dt_save end-1/dt_save:end]) = [];
Zfilt(:, [1:1/dt_save end-1/dt_save:end]) = [];

% Define thresholds as 5*STD of the power in each band in each area
sigma_std = 5;
deltathr = sigma_std*std(Zdelta, [], 2);
thetathr = sigma_std*std(Ztheta, [], 2);
alphathr = sigma_std*std(Zalpha, [], 2);
betathr = sigma_std*std(Zbeta, [], 2);
filtthr = sigma_std*std(Zfilt,[],2);
clear Zalpha Zbeta Ztheta Zfilt


%% 2 - Load data in the working point and filter into bands

% Change the "simu" to look at other scenarios outside the optimal 
% working point 

load(fullfile(folder,simulation_file{simu}));
N = size(Zsave, 1);

Zsave = Zsave./(5*std(Zsave(:)));
Zsave = Zsave(Order, :);
Zdelta = zeros(size(Zsave));
Ztheta = zeros(size(Zsave));
Zalpha = zeros(size(Zsave));
Zbeta = zeros(size(Zsave));
Zfilt  = zeros(size(Zsave));

% Bandpass filter (via FFT) to obtain power at each frequency band 
% Cycle through nodes 
for n = 1 : N
    
    % Syntax is bandpasshopf( data, [f_min f_max], fs ) 
    Zdelta(n,:) = bandpasshopf(Zsave(n, :), delta, 1/dt_save);
    Ztheta(n,:) = bandpasshopf(Zsave(n, :), theta, 1/dt_save);
    Zalpha(n,:) = bandpasshopf(Zsave(n, :), alpha, 1/dt_save);
    Zbeta(n,:)  = bandpasshopf(Zsave(n, :), beta, 1/dt_save);
    Zfilt(n,:)  = bandpasshopf(Zsave(n, :), [0.01 low_pass], 1/dt_save);
    
end % nodes 
clear Zsave

% Remove the first and last second after bandpass filtering
Zdelta(:, [1:1/dt_save end-1/dt_save:end]) = [];
Ztheta(:, [1:1/dt_save end-1/dt_save:end]) = [];
Zalpha(:, [1:1/dt_save end-1/dt_save:end]) = [];
Zbeta(:, [1:1/dt_save end-1/dt_save:end]) = [];
Zfilt(:, [1:1/dt_save end-1/dt_save:end]) = [];

% Calculate amplitude envelopes for each band,
% using the Hilbert transform 
Env_Delta = abs(hilbert(Zdelta'))';
Env_Theta = abs(hilbert(Ztheta'))';
Env_Alpha = abs(hilbert(Zalpha'))';
Env_Beta  = abs(hilbert(Zbeta'))';
Env_Zfilt = abs(hilbert(Zfilt'))';

clear Zdelta Zalpha Zbeta Ztheta

% Detect formation of a MOMs in each band
% Criteria: a node engages in a MOM if the amplitude increases
% 5 standard deviations (threshold) above the baseline amplitude 
% in that frequency band/range
T_Delta = Env_Delta > repmat(deltathr, 1, size(Env_Delta, 2));
T_Theta = Env_Theta > repmat(thetathr, 1, size(Env_Theta, 2));
T_Alpha = Env_Alpha > repmat(alphathr, 1, size(Env_Alpha, 2));
T_Beta = Env_Beta > repmat(betathr, 1, size(Env_Beta, 2)); 

% clear Env_Delta Env_Theta Env_Alpha Env_Beta


%% 3 - Calculate Shannon entropy of covariance eigenvalues 
%
% Calculate the Shannon entropy of the eigenvalue distribution of the 
% phase covariance of the filtered data throughout time (sliding window) 

% Define window and step size 
% of the sliding window 
win_size = 100;
win_step = 1;

% Compute the instantaneous phase of the filtered signal 
Phase_filt = angle(hilbert(Zfilt'))';

% Compute the order parameter of the filtered signal 
OP = abs(mean(exp(1i*Phase_filt)));
OP_mean = mean(OP);
OP = bandpasshopf(OP-OP_mean,[0.0001 3], 1/dt_save) + OP_mean; 

% Initialize variables to store entropy and probability distributions
total_time_points = size(Zfilt, 2);
n_wins = floor((total_time_points - win_size) / win_step);
Shan_Entropy = zeros(1, n_wins); 
Eigen_Prob_Distrib = zeros(size(Zfilt, 1), n_wins); 

% Calculate the Shannon entropy for each sliding window
for idx = 1 : n_wins 
    center_point = (idx - 1) * win_step + floor(win_size / 2) + 1;
    win_start = center_point - floor(win_size / 2);
    win_end = center_point + floor(win_size / 2);

    win_inds = win_start : win_end;

    % Covariance and eigenvalue distribution
    cov_matrix = cov(Phase_filt(:, win_inds)');
    eigenvalues = eig(cov_matrix);
    eigenvalues_norm = eigenvalues / sum(eigenvalues);

    % Store probability distribution and calculate entropy
    Eigen_Prob_Distrib(:, idx) = eigenvalues_norm;
    Shan_Entropy(idx) = -sum(eigenvalues_norm .* log2(eigenvalues_norm));
end

%% 4 - Calculate standard deviation of SE over time 
%

SE_std = std(Shan_Entropy);
save(fullfile(folder_pars, 'SE_std.mat'), "SE_std", '-mat');

%% 5 - Calculate temporal correlation between MAE and SE
%
% Calculate the temporal correlation between the Mean Amplitude Envelope 
% and the Shannon Entropy + corresponding statistics 

Mean_Env_Zfilt = mean(Env_Zfilt); 

% Compute and save correlation 
[SE_MAE_corr.cc, SE_MAE_corr.p] = corrcoef(Shan_Entropy, ...
    Mean_Env_Zfilt((win_size/2) + 1 : end - (win_size/2)));
save(fullfile(folder_pars, 'SE_MAE_corr.mat'), "SE_MAE_corr", '-mat');


%% 6 - Plot of signals over time

Time_to_plot = 10; % Plot only first 25s

load AAL_labels.mat label90

fig = figure('Color', 'white');
fig.Position(3:4) = fig.Position(3:4)*8;
hold on

% Delta patches
for n=1:N
    u=0;
    y=[];    x=[];
    for tx=find(T_Delta(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,[248 247 85]/255);
    %set(p,'LineStyle','none','FaceColor','y','FaceAlpha',0.3);
    set(p,'LineStyle','none','FaceColor',[248 247 85]/255,'FaceAlpha',0.6);
end
ylim([0 N+1])

% Theta patches
for n=1:N
    u=0;
    y=[];    x=[];
    for tx=find(T_Theta(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,[238 188 86]/255);
    %set(p,'LineStyle','none','FaceColor','r','FaceAlpha',0.3);
    set(p,'LineStyle','none','FaceColor',[238 188 86]/255,'FaceAlpha',0.6);
end
ylim([0 N+1])

% Alpha patches
for n=1:N
    u=0;
    y=[];    x=[];
    for tx=find(T_Alpha(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,[77 153 226]/255);
    %set(p,'LineStyle','none','FaceColor','b','FaceAlpha',0.3);
    set(p,'LineStyle','none','FaceColor',[77 153 226]/255,'FaceAlpha',0.6);    
end

% Beta patches
for n=1:N
    u=0;
    y=[];
    x=[];
    for tx=find(T_Beta(n,1:Time_to_plot/dt_save))
        u=u+1;
        y(:,u)= [n-0.5 n-0.5  n+.5 n+.5];
        x(:,u) = [tx-1 tx tx tx-1];
    end
    p=patch(x.*dt_save,y,[131 203 117]/255);
    %set(p,'LineStyle','none','FaceColor','g','FaceAlpha',0.3);
    set(p,'LineStyle','none','FaceColor','#83CB75','FaceAlpha',0.6);  
end

plot(0:dt_save:(length(Zfilt)-1)*dt_save,(1:N)'.*ones(size(Zfilt))+(Zfilt./(2*filtthr)),'k')
xlabel('Time (seconds)', 'FontSize', 18, 'FontName', 'Helvetica')
box off
set(gca,'YTick',1:N,'Fontsize',16)
set(gca,'YTickLabel',[])
set(gca,'YTickLabel',label90(Order,:))
ylim([0 N+1])
xlim([0 Time_to_plot])
saveas(gcf, fullfile(folder_img, 'MOM'), 'png');

% Mean Amplitude Envelope Correlation, OP, Entropy  
OP = OP(1 : (Time_to_plot/dt_save) + 1);
Mean_Env_Zfilt = Mean_Env_Zfilt(:, 1 : (Time_to_plot/dt_save) + 1);
Shan_Entropy = Shan_Entropy(1 : length(Mean_Env_Zfilt) - (win_size/2));
Eigen_Prob_Distrib = Eigen_Prob_Distrib(:, 1 : length(Mean_Env_Zfilt) - ((win_size/2) + 1 ));
Time_Mean_Env_Zfilt = 0: dt_save : (length(Mean_Env_Zfilt) - 1)*dt_save;

% Mean Amplitude Envelope and Kuramoto Order Parameter
fig = figure();
fig.Position(3:4) = fig.Position(3:4)*8;
[AX,H1,H2] = plotyy(0 : dt_save : (length(OP)-1)*dt_save, Mean_Env_Zfilt, ...
    0 : dt_save : (length(OP)-1)*dt_save, OP);
set(AX(2), 'ytick', 0:0.5:1, 'ylim', [0 1], 'xlim', [0 Time_to_plot], 'Fontsize', 16, 'YColor', [238 188 86]/255)
set(AX(1), 'ytick', [0.1 0.2], 'ylim', [0.05 0.25], 'xlim', [0 Time_to_plot], 'Fontsize', 16, 'YColor', [77 153 226]/255)
set(H1, 'Color', [77 153 226]/255, 'LineWidth', 2); 
set(H2, 'Color', [238 188 86]/255, 'LineWidth', 2); 
xlabel('Time (seconds)', 'Fontsize', 16)
legend({'Mean Amplitude Envelope', 'Kuramoto Order Parameter'}, ...
    'Orientation', 'horizontal', 'Fontsize', 16)
box off
saveas(gcf, fullfile(folder_img, 'mean_AE_OP'), 'png');

% Mean Amplitude Envelope and Shannon Entropy  
fig = figure();
fig.Position(3:4) = fig.Position(3:4)*8;
[AX,H1,H2] = plotyy(Time_Mean_Env_Zfilt((win_size/2) + 1 : end), ...
    Mean_Env_Zfilt((win_size/2) + 1 : end), ...
    Time_Mean_Env_Zfilt((win_size/2) + 1 : end), Shan_Entropy);
set(AX(2), 'ytick', 3:5, 'ylim', [3 5], 'xlim', [0 Time_to_plot], 'Fontsize', 16, 'YColor', [238 188 86]/255)
set(AX(1), 'ytick', [0.1 0.2], 'ylim', [0.05 0.25], 'xlim', [0 Time_to_plot], 'Fontsize', 16, 'YColor', [77 153 226]/255)
set(H1, 'Color', [77 153 226]/255, 'LineWidth', 2); 
set(H2, 'Color', [238 188 86]/255, 'LineWidth', 2); 
xlabel('Time (seconds)', 'Fontsize', 16)
legend({'Mean Amplitude Envelope', 'Shannon Entropy'}, ...
   'Orientation', 'horizontal', 'Fontsize', 16)
box off
saveas(gcf, fullfile(folder_img, 'mean_AE_SE'), 'png');

% Eigenvector Probability Distribution 
fig = figure();
fig.Position(3:4) = fig.Position(3:4)*8;
im = imagesc(Time_Mean_Env_Zfilt((win_size/2) + 1 : end), ...
    1 : size(Eigen_Prob_Distrib, 1), (Eigen_Prob_Distrib.^(1/2)));
xlabel('Time (seconds)', 'Fontsize', 16); im.Interpolation = 'bilinear';
c = colorbar; c.TickLabels = '';
%ylabel(c, 'Eigenvalue','FontSize', 12, 'Rotation',270);
saveas(gcf, fullfile(folder_img, 'Eigen_Prob_Distrib'), 'png');

% Barplots of minimum and maximum mean Eigen. Probab. Distribution 
[mini, idx_min] = min(Eigen_Prob_Distrib(end, :));
[maxi, idx_max] = max(Eigen_Prob_Distrib(end, :));

figure; % Minimum 
bar(Eigen_Prob_Distrib(:, idx_min), 'FaceColor', [77 153 226]/255);
saveas(gcf, fullfile(folder_img, 'Eigen_Prob_Distrib_min'), 'png');

figure; % Maximum 
bar(Eigen_Prob_Distrib(:, idx_max), 'FaceColor', [77 153 226]/255);
saveas(gcf, fullfile(folder_img, 'Eigen_Prob_Distrib_max'), 'png');
