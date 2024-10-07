%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to detect MOMs (Metastable Oscillatory Modes) in different 
% frequency bands for a selected point in the parameter space and calculate
% Shannon Entropy of the LFP phase covariance matrices throughout time 
%
% Needs: bandpasshopf.m, conver_back_to_time.m, subplot_tight.m
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
%% Initialization and Parameters
simu = 2; % Select your simulation of interest
simulation_names = {'NoDelays', 'WeakK', 'IntermediateK', 'StrongK', 'LongDelays'};

% Frequency bands to analyze (in Hz)
delta = [0.5 4];
theta = [4 8];
alpha = [8 13];
beta = [13 30];
low_pass = 30;

folder = 'C:\Users\marta\OneDrive\Documentos\LASEEB\Simulations\FEP_MOMs\';

K = '1E1'; MD = '0';
simulation_file = {strcat('Hopf_Simu_K', K, '_MD0.mat'), strcat('Hopf_Simu_K', K, '_MD', MD, '.mat')};

folder_pars = fullfile(folder, strcat('K', K, '_MD', MD)); mkdir(folder_pars);
folder_img = fullfile(folder_pars, 'Imgs'); mkdir(folder_img);

%% 1 - Load Baseline Data and Define Thresholds
% Load data for baseline (MD = 0)
% Define thresholds based on 5*STD of the power in each band in each area

load(fullfile(folder, simulation_file{1})); % Baseline with MD = 0

% Variables 
N = size(Zsave, 1);
Order = [1:2:N N:-2:2];
Zsave = Zsave./(5*std(Zsave(:)));
Zsave = Zsave(Order, :);
Zbands = {'delta', 'theta', 'alpha', 'beta', 'filt'};
Zfiltered = struct();

% Apply bandpass filter to obtain power in each frequency band
for n = 1:N
    Zfiltered.delta(n, :) = bandpasshopf(Zsave(n, :), delta, 1/dt_save);
    Zfiltered.theta(n, :) = bandpasshopf(Zsave(n, :), theta, 1/dt_save);
    Zfiltered.alpha(n, :) = bandpasshopf(Zsave(n, :), alpha, 1/dt_save);
    Zfiltered.beta(n, :) = bandpasshopf(Zsave(n, :), beta, 1/dt_save);
    Zfiltered.filt(n, :) = bandpasshopf(Zsave(n, :), [0.01 low_pass], 1/dt_save);
end
clear Zsave

% Remove the first and last second after bandpass filtering
fields = fieldnames(Zfiltered);
for i = 1:numel(fields)
    Zfiltered.(fields{i})(:, [1:1/dt_save end-1/dt_save:end]) = [];
end

% Define thresholds for each band and region 
sigma_std = 5;
thresholds = struct();
for i = 1:numel(fields)
    thresholds.(fields{i}) = sigma_std * std(Zfiltered.(fields{i}), [], 2);
end
clear Zfiltered

%% 2 - Load Working Point Data and Filter into Bands
% Load data at the working point and filter into frequency bands

load(fullfile(folder, simulation_file{simu}));

% Variables 
N = size(Zsave, 1);
Zsave = Zsave./(5*std(Zsave(:)));
Zsave = Zsave(Order, :);

% Apply bandpass filter to obtain power in each frequency band
for n = 1:N
    Zfiltered.delta(n, :) = bandpasshopf(Zsave(n, :), delta, 1/dt_save);
    Zfiltered.theta(n, :) = bandpasshopf(Zsave(n, :), theta, 1/dt_save);
    Zfiltered.alpha(n, :) = bandpasshopf(Zsave(n, :), alpha, 1/dt_save);
    Zfiltered.beta(n, :) = bandpasshopf(Zsave(n, :), beta, 1/dt_save);
    Zfiltered.filt(n, :) = bandpasshopf(Zsave(n, :), [0.01 low_pass], 1/dt_save);
end
clear Zsave

% Remove the first and last second after bandpass filtering
for i = 1:numel(fields)
    Zfiltered.(fields{i})(:, [1:1/dt_save end-1/dt_save:end]) = [];
end

% Calculate amplitude envelopes for each band using the Hilbert transform
Env = struct();
for i = 1:numel(fields) 
    Env.(fields{i}) = abs(hilbert(Zfiltered.(fields{i})'))';
end

% Detect formation of MOMs in each band
T = struct();
for i = 1:numel(fields) - 1
    T.(fields{i}) = Env.(fields{i}) > repmat(thresholds.(fields{i}), 1, size(Env.(fields{i}), 2));
end

%% 3 - Calculate MOM Durations, Size, Coalition Size, and Occupancy
% Calculate MOM durations and size, coalition size, and occupancy for 
% each frequency band

% Initialize structures for storing results
MOM_Durations = struct();
MOM_Mean_Duration = struct();
MOM_Std_Duration = struct();

Coalition_Members = struct();
Coalition_Members_Mean = struct();
Coalition_Members_Std = struct();

MOM_Occupancy = struct();

% Loop over each frequency band to calculate MOM metrics
for i = 1:4
    band_name = fields{i};
    
    durations = [];
    for n = 1:N
        % Detect switches in and out of this state
        a = find(diff(T.(band_name)(n, :)) == 1);  % in (0-->1)
        b = find(diff(T.(band_name)(n, :)) == -1); % out (1-->0)
        
        % Discard the cases where state starts or ends 'in'
        if length(b) > length(a)
            b(1) = [];
        elseif length(a) > length(b)
            a(end) = [];
        elseif ~isempty(a) && ~isempty(b) && a(1) > b(1)
            b(1) = [];
            a(end) = [];
        end

        if ~isempty(a) && ~isempty(b)
            durations = [durations b-a];
        end        
    end

    % MOM durations 
    MOM_Durations.(band_name) = durations * dt_save;
    MOM_Mean_Duration.(band_name) = mean(durations) * dt_save;
    MOM_Std_Duration.(band_name) = std(durations) * dt_save;
    
    % Coalition size
    members = sum(T.(band_name)(:, 2:end-2));
    Coalition_Members.(band_name) = members(members > 0);
    Coalition_Members_Mean.(band_name) = mean(Coalition_Members.(band_name));
    Coalition_Members_Std.(band_name) = std(Coalition_Members.(band_name));
    
    % MOM occupancy
    MOM_Occupancy.(band_name) = sum(sum(T.(band_name)(:, 2:end-2))) / numel(T.(band_name));
end

% Define the threshold for the number of simultaneous coalition members 
% in order to detect a MOM 
MOM_Occurences = struct();
threshold_members = struct();   % since the distribution of coalition members 
                                % (for each frequency band) approximates an
                                % exponential distribution, we can fit an 
                                % exponential equation to the distribution 
                                % and find its tail

% Find one threshold for all bands 
data = [Coalition_Members.delta Coalition_Members.theta ...
    Coalition_Members.alpha Coalition_Members.beta];
lambda_hat = 1 / mean(data); p = 0.5;
threshold_members_allbands = round(-log(1 - p) / lambda_hat); %50th percentile 

for i = 1:4
    band_name = fields{i};
    
    data = Coalition_Members.(band_name);
    
    % Fit the exponential distribution
    % Estimate the parameter lambda (rate parameter)
    lambda_hat = 1 / mean(data); % MLE for the exponential distribution
    
    % Determine the tail threshold
    % Define the desired percentile
    p = 0.90;
    
    % Calculate the threshold T for the 90th percentile
    %threshold_members.(band_name) = round(-log(1 - p) / lambda_hat);   % Option 1
    threshold_members.(band_name) = threshold_members_allbands;         % Option 2

    % Apply threshold 
    % Identify time points where the number of engaged nodes exceeds the threshold
    MOM_times = sum(T.(band_name)) > threshold_members.(band_name);
    
    % Store the results
    MOM_Occurences.(band_name) = MOM_times;

end

%% 4 - Plot MOM Size Distribution 
% Plot MOM size distribution for each frequency band

if sum(T.(band_name))>=1
    % Plot MOM size distribution (without size threshold)
    figure('Position', [100, 100, 1200, 600]);
    x_label = 'MOM Size';
    y_label = 'Frequency';
    
    bands = {'Delta', 'Theta', 'Alpha', 'Beta'};
    for i = 1:4
        band_name = fields{i};
    
        MOM_sizes = sum(T.(band_name), 1);
        x_ticks = min(MOM_sizes):1:max(MOM_sizes);
        generate_labels = @(ticks) arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
        x_tick_labels = generate_labels(x_ticks);
        x_tick_labels(2:2:end) = {''}; 
    
        subplot(2, 2, i); 
        histogram(MOM_sizes, 'BinEdges', x_ticks - 0.5); hold on;
        
        % Plot vertical line with the threshold 
        xline(threshold_members.(band_name) - 0.5, 'r', 'LineWidth', 2); 
    
        title([bands{i}, ' Coalition Size Distribution']); 
        xlabel(x_label); ylabel(y_label);
        xticks(x_ticks); xtickangle(45); 
        xticklabels(x_tick_labels);
        set(gca, 'FontSize', 8); xlim([(min(x_ticks) - 0.5) max(x_ticks)]);
        grid on;
    end
    saveas(gcf, fullfile(folder_pars, 'MOM_size_distribution_untrhesh'), 'png');
end

if sum(T.(band_name).*MOM_Occurences.(band_name))>=1
    % Plot MOM size distribution (w/ size threshold)
    figure('Position', [100, 100, 1200, 600]);
    x_label = 'MOM Size';
    y_label = 'Frequency';
    
    bands = {'Delta', 'Theta', 'Alpha', 'Beta'};
    for i = 1:4
        band_name = fields{i};
    
        MOM_sizes = sum(T.(band_name).*MOM_Occurences.(band_name));
        MOM_sizes(MOM_sizes==0)=[];
    
        x_ticks = min(MOM_sizes):1:max(MOM_sizes);
        generate_labels = @(ticks) arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
        x_tick_labels = generate_labels(x_ticks);
        x_tick_labels(2:2:end) = {''}; 
    
        subplot(2, 2, i); 
        histogram(MOM_sizes, 'BinEdges', x_ticks - 0.5);
        title([bands{i}, ' Coalition Size Distribution']); 
        xlabel(x_label); ylabel(y_label);
        xticks(x_ticks); xtickangle(45); 
        xticklabels(x_tick_labels);
        set(gca, 'FontSize', 8); xlim([(min(x_ticks) - 0.5) max(x_ticks)]);
        grid on;
    end
    saveas(gcf, fullfile(folder_pars, 'MOM_size_distribution'), 'png');
end

%% 5 - Calculate Shannon Entropy of Covariance Eigenvalues, OP, Metastability
% Calculate the Shannon entropy of the eigenvalue distribution of the 
% phase covariance of the filtered data throughout time (sliding window)

% Define window and step size
win_size = 100;
win_step = round(0.25*win_size);

% Compute the instantaneous phase of the filtered signal 
Phase_filt = angle(hilbert(Zfiltered.filt'))';

% Initialize variables to store entropy and probability distributions
total_time_points = size(Zfiltered.filt, 2);
n_wins = floor((total_time_points - win_size) / win_step);
Shan_Entropy = zeros(1, n_wins); 
Eigen_Prob_Distrib = zeros(size(Zfiltered.filt, 1), n_wins); 

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

% Compute the Order Parameter (OP)
OP = abs(mean(exp(1i * Phase_filt)));
OP_mean = mean(OP);
OP = bandpasshopf(OP - OP_mean, [0.0001 3], 1/dt_save) + OP_mean; 

SE_std = std(Shan_Entropy);               % Standard deviation of SE
Metastability = std(OP);                  % Metastability (standard 
                                          % deviation of OP) 

% Save variables 
save(fullfile(folder_pars, 'Shan_Entropy.mat'), "Shan_Entropy", '-mat');
save(fullfile(folder_pars, 'SE_std.mat'), "SE_std", '-mat');
save(fullfile(folder_pars, 'Metastability.mat'), 'Metastability', '-mat');

%% 6 - Calculate Temporal Correlation Between Mean Amplitude Envelope and Shannon Entropy
% Average power in windows before computing correlation 
Env_Zfilt_win = zeros(size(Zfiltered.filt, 1), n_wins);
for idx = 1 : n_wins 
    center_point = (idx - 1) * win_step + floor(win_size / 2) + 1;
    win_start = center_point - floor(win_size / 2);
    win_end = center_point + floor(win_size / 2);

    win_inds = win_start : win_end;
    Env_Zfilt_win(:, idx) = mean(Env.filt(:, win_inds), 2);
end

% Compute and save correlation 
Mean_Env_Zfilt = mean(Env_Zfilt_win);
[SE_MAE_corr.cc, SE_MAE_corr.p] = corrcoef(Shan_Entropy, Mean_Env_Zfilt);
save(fullfile(folder_pars, 'SE_MAE_corr.mat'), "SE_MAE_corr", '-mat');

%% 7 - Calculate Temporal Correlation Between Total Coalition Size and Shannon Entropy
% Calculate the total coalition size across all frequency bands and correlate
% it with the Shannon entropy.

% Initialize variable to store total coalition size across all bands
Total_Coalition_Size = zeros(1, size(T.delta, 2));

% Sum the coalition sizes across all frequency bands at each time point
for i = 1:4
    band_name = fields{i};
    
    % Add the coalition size for each band to the total coalition size
    Total_Coalition_Size = Total_Coalition_Size + sum(T.(band_name));
end

% Windowed 
Total_Coalition_Size_win = zeros(size(Total_Coalition_Size, 1), n_wins);
for idx = 1 : n_wins 
    center_point = (idx - 1) * win_step + floor(win_size / 2) + 1;
    win_start = center_point - floor(win_size / 2);
    win_end = center_point + floor(win_size / 2);

    win_inds = win_start : win_end;
    Total_Coalition_Size_win(:, idx) = mean(Total_Coalition_Size(:, win_inds), 2);
end
Total_Coalition_Size = Total_Coalition_Size_win;

% Compute correlation with Shannon Entropy
[SE_Coalition_Size_corr.cc, SE_Coalition_Size_corr.p] = corrcoef(Shan_Entropy, Total_Coalition_Size);

% Store results
save(fullfile(folder_pars, 'SE_Coalition_Size_corr.mat'), 'SE_Coalition_Size_corr', '-mat');

%% 8 - Calculate Temporal Correlation Between Mean Amplitude Envelope and Shannon Entropy 
% Pool all time points in which MOMs occur across any frequency band and 
% those in which MOMs don't occur.

% Identify MOM and no-MOM time points across all bands
% All_MOM_time_points = false(1, size(T.delta, 2));
% 
% % Loop over each frequency band to accumulate MOM occurrences
% for i = 1:4
%     band_name = fields{i};
%     All_MOM_time_points = All_MOM_time_points | MOM_Occurences.(band_name);
% end
% 
% % Define time points for MOMs and no-MOMs
% MOM_time_points = All_MOM_time_points(win_size/2 : end - win_size/2 - 1);
% no_MOM_time_points = ~All_MOM_time_points(win_size/2 : end - win_size/2 - 1);
% 
% % Mean amplitude envelope for MOMs
% Mean_Env_MOMs = Mean_Env_Zfilt(MOM_time_points);
% Shan_Entropy_MOMs = Shan_Entropy(MOM_time_points);
% 
% % Mean amplitude envelope for no-MOMs
% Mean_Env_noMOMs = Mean_Env_Zfilt(no_MOM_time_points);
% Shan_Entropy_noMOMs = Shan_Entropy(no_MOM_time_points);
% 
% % Initialize structures to store correlation results
% SE_MAE_corr_MOMs = struct();
% SE_MAE_corr_noMOMs = struct();
% 
% % Compute correlation for MOMs
% if ~isempty(Mean_Env_MOMs) && length(Mean_Env_MOMs) > 1
%     [SE_MAE_corr_MOMs.cc, SE_MAE_corr_MOMs.p] = corrcoef(Shan_Entropy_MOMs, Mean_Env_MOMs);
% else
%     SE_MAE_corr_MOMs.cc = NaN;
%     SE_MAE_corr_MOMs.p = NaN;
% end
% 
% % Compute correlation for no-MOMs
% if ~isempty(Mean_Env_noMOMs) && length(Mean_Env_noMOMs) > 1
%     [SE_MAE_corr_noMOMs.cc, SE_MAE_corr_noMOMs.p] = corrcoef(Shan_Entropy_noMOMs, Mean_Env_noMOMs);
% else
%     SE_MAE_corr_noMOMs.cc = NaN;
%     SE_MAE_corr_noMOMs.p = NaN;
% end
% 
% % Store results
% save(fullfile(folder_pars, 'SE_MAE_corr_MOMs.mat'), 'SE_MAE_corr_MOMs', '-mat');
% save(fullfile(folder_pars, 'SE_MAE_corr_noMOMs.mat'), 'SE_MAE_corr_noMOMs', '-mat');

%% 9.1 - Plot Power (<30Hz) time-series and MOMs (without size threshold)
% Power (<30Hz) = Mean Amplitude Envelope (<30Hz)

Time_to_plot = 10; % Plot only first 10s

load AAL_labels.mat label90

fig = figure('Color', 'white');
fig.Position(3:4) = fig.Position(3:4)*8;
hold on

% Define patch colors 
patch_color.delta = [248 247 85]/255;
patch_color.theta = [238 188 86]/255;
patch_color.alpha = [77 153 226]/255;
patch_color.beta = [131 203 117]/255;

for i = 1 : 4
    band_name = fields{i};
    for n=1:N
        u=0;
        y=[];    x=[];
        for tx=find(T.(band_name)(n, 1:Time_to_plot/dt_save))
            u = u + 1;
            y(:, u)= [n-0.5 n-0.5  n+.5 n+.5];
            x(:, u) = [tx-1 tx tx tx-1];
        end
        p=patch(x.*dt_save, y, patch_color.(band_name));
        set(p, 'LineStyle', 'none', 'FaceColor', patch_color.(band_name), 'FaceAlpha', 0.6);
    end
    ylim([0 N+1])
end

plot(0:dt_save:(length(Zfiltered.filt)-1)*dt_save, ...
    (1:N)'.*ones(size(Zfiltered.filt))+(Zfiltered.filt./(2*thresholds.filt)), 'k')
xlabel('Time (seconds)', 'FontSize', 18, 'FontName', 'Helvetica')
box off
set(gca, 'YTick', 1:N, 'Fontsize', 16)
set(gca, 'YTickLabel', [])
set(gca, 'YTickLabel', label90(Order, :))
ylim([0 N+1])
xlim([0 Time_to_plot])
saveas(gcf, fullfile(folder_img, 'MOM_untresh'), 'png');


%% 9.1 - Plot Power (<30Hz) time-series and MOMs (w/ size threshold)
% Power (<30Hz) = Mean Amplitude Envelope (<30Hz)

fig = figure('Color', 'white');
fig.Position(3:4) = fig.Position(3:4)*8;
hold on

% Define patch colors 
patch_color.delta = [248 247 85]/255;
patch_color.theta = [238 188 86]/255;
patch_color.alpha = [77 153 226]/255;
patch_color.beta = [131 203 117]/255;

for i = 1 : 4
    band_name = fields{i};
    for n=1:N
        u=0;
        y=[];    x=[];
        data = T.(band_name).*MOM_Occurences.(band_name);
        for tx=find(data(n, 1:Time_to_plot/dt_save))
            u = u + 1;
            y(:, u)= [n-0.5 n-0.5  n+.5 n+.5];
            x(:, u) = [tx-1 tx tx tx-1];
        end
        p=patch(x.*dt_save, y, patch_color.(band_name));
        set(p, 'LineStyle', 'none', 'FaceColor', patch_color.(band_name), 'FaceAlpha', 0.6);
    end
    ylim([0 N+1])
end

plot(0:dt_save:(length(Zfiltered.filt)-1)*dt_save, ...
    (1:N)'.*ones(size(Zfiltered.filt))+(Zfiltered.filt./(2*thresholds.filt)), 'k')
xlabel('Time (seconds)', 'FontSize', 18, 'FontName', 'Helvetica')
box off
set(gca, 'YTick', 1:N, 'Fontsize', 16)
set(gca, 'YTickLabel', [])
set(gca, 'YTickLabel', label90(Order, :))
ylim([0 N+1])
xlim([0 Time_to_plot])
saveas(gcf, fullfile(folder_img, 'MOM'), 'png');

%% 10 - Plot Power(<30Hz), Shannon Entropy and Coalition Size time-series 
% Power (<30Hz) = Mean Amplitude Envelope (<30Hz)

% Adjust
SElim = [2 6];
MAElim = [0.01 0.3];
CoSizelim = [0 210];

SElim = [floor(min(Shan_Entropy)) ceil(max(Shan_Entropy))];
CoSizelim = [floor(min(Total_Coalition_Size)) ceil(max(Total_Coalition_Size))];

% Align time-series to plot 
%Mean_Env_Zfilt = Mean_Env_Zfilt(1 : (Time_to_plot/(dt_save*win_step)) + 1);
%Shan_Entropy = Shan_Entropy(1 : (Time_to_plot/(dt_save*win_step)) + 1);
%Total_Coalition_Size = Total_Coalition_Size(1 : (Time_to_plot/(dt_save*win_step)) + 1);
%Eigen_Prob_Distrib = Eigen_Prob_Distrib(:, 1 : (Time_to_plot/(dt_save*win_step)) + 1);
Time_Mean_Env_Zfilt = 0: (dt_save*win_step) : (length(Mean_Env_Zfilt) - 1)*(dt_save*win_step);

% Mean Amplitude Envelope and Shannon Entropy  
fig = figure();
fig.Position(3:4) = fig.Position(3:4)*6;
[AX,H1,H2] = plotyy(Time_Mean_Env_Zfilt, Mean_Env_Zfilt, Time_Mean_Env_Zfilt, Shan_Entropy);
set(AX(1), 'ytick', MAElim, 'ylim', MAElim, 'xlim', [0 Time_to_plot], 'Fontsize', 36, 'YColor', [77 153 226]/255);
set(AX(2), 'ytick', SElim, 'ylim', SElim, 'xlim', [0 Time_to_plot], 'Fontsize', 36, 'YColor', [238 188 86]/255);
AX(1).Position = AX(1).Position + [0, 0.05, 0, -0.05]; 
AX(2).Position = AX(2).Position + [0, 0.05, 0, -0.05]; 
set(H1, 'Color', [77 153 226]/255, 'LineWidth', 3.5); 
set(H2, 'Color', [238 188 86]/255, 'LineWidth', 3.5); 
%legend('Power 30 Hz', 'Shannon Entropy');
xlabel('Time (seconds)', 'Fontsize', 36);
box off
saveas(gcf, fullfile(folder_img, 'mean_MAE_SE'), 'png');

% Coalition Size and Shannon Entropy
fig = figure();
fig.Position(3:4) = fig.Position(3:4)*6;
[AX,H1,H2] = plotyy(Time_Mean_Env_Zfilt, Total_Coalition_Size, Time_Mean_Env_Zfilt, Shan_Entropy);
set(AX(1), 'ytick', CoSizelim, 'ylim', CoSizelim, 'xlim', [0 47.7], 'Fontsize', 36, 'YColor', [77 153 226]/255)
set(AX(2), 'ytick', SElim, 'ylim', SElim, 'xlim', [0 47.7], 'Fontsize', 36, 'YColor', [238 188 86]/255)
AX(1).Position = AX(1).Position + [0, 0.05, 0, -0.05]; 
AX(2).Position = AX(2).Position + [0, 0.05, 0, -0.05]; 
set(H1, 'Color', [77 153 226]/255, 'LineWidth', 2.5); 
set(H2, 'Color', [238 188 86]/255, 'LineWidth', 2.5); 
%legend('Coalitions Size', 'Shannon Entropy', 'Orientation', 'Vertical');
xlabel('Time (seconds)', 'Fontsize', 36);
box off
saveas(gcf, fullfile(folder_img, 'CoSize_SE'), 'png');

% Eigenvector Probability Distribution 
fig = figure();
fig.Position(3:4) = fig.Position(3:4)*8;
im = imagesc(Time_Mean_Env_Zfilt, 1 : size(Eigen_Prob_Distrib, 1), (Eigen_Prob_Distrib.^(1/2)));
xlabel('Time (seconds)', 'Fontsize', 16); im.Interpolation = 'bilinear';
c = colorbar; c.TickLabels = '';
%ylabel(c, 'Eigenvalue', 'FontSize', 12, 'Rotation', 270);
saveas(gcf, fullfile(folder_img, 'Eigen_Prob_Distrib'), 'png');

% Barplots of minimum and maximum mean Eigen. Probab. Distribution 
[~, idx_min] = min(Eigen_Prob_Distrib(end, :));
[~, idx_max] = max(Eigen_Prob_Distrib(end, :));

figure; % Minimum 
bar(Eigen_Prob_Distrib(:, idx_min), 'FaceColor', [77 153 226]/255);
saveas(gcf, fullfile(folder_img, 'Eigen_Prob_Distrib_min'), 'png');

figure; % Maximum 
bar(Eigen_Prob_Distrib(:, idx_max), 'FaceColor', [77 153 226]/255);
saveas(gcf, fullfile(folder_img, 'Eigen_Prob_Distrib_max'), 'png');
