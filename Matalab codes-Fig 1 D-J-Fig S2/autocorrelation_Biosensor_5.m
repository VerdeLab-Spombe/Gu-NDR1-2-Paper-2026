%% Temporal autocorrelation from CSV (rows 2:26 = samples, cols 2:100 = time)
clear; clc;

% ---------------- user params ----------------
infile = 'T-Windows_Ave_Example_1.csv';  % your input CSV
dt_sec = 10;       % sampling interval
lag_sec = 200;      % total lag to each side (±200 s)
% ---------------------------------------------

% Read CSV and extract the exact numeric block
M = readmatrix(infile);            % reads numeric content from CSV
X = M(2:50, 2:101);                % rows=samples, cols=time
[numSamples, numTime] = size(X);   %#ok<NASGU>

% Clean & zero-mean each sample (row)
for i = 1:numSamples
    xi = X(i, :).';
    if any(isnan(xi))
        xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
    end
    xi = detrend(xi, 'constant');  % remove mean
    X(i, :) = xi.';
end

% Lags (±200 s at 10 s step => ±20 samples)
maxlag = round(lag_sec / dt_sec);      % 20
lag_samples_vec = (-maxlag:maxlag).';
lag_seconds_vec = lag_samples_vec * dt_sec;
nLags = numel(lag_seconds_vec);

% ACF per sample (row-wise)
ACF = zeros(nLags, numSamples);
for i = 1:numSamples
    [acfi, ~] = xcorr(X(i, :).', maxlag, 'coeff');
    ACF(:, i) = acfi(:);
end

% Mean ± SEM across samples
ACF_mean = mean(ACF, 2, 'omitnan');
ACF_sem  = std(ACF, 0, 2, 'omitnan') / sqrt(numSamples);

% ---------- Write CSV outputs ----------
% Lags
T_lags = table(lag_seconds_vec, lag_samples_vec, ...
    'VariableNames', {'lag_seconds','lag_samples'});
writetable(T_lags, 'lags_T-Windows_Ave_Example_1.csv');

% Per-sample ACF (columns: Sample_1 ... Sample_25)
sampleNames = "Sample_" + string(1:numSamples);
T_per = array2table([lag_seconds_vec, ACF], ...
    'VariableNames', [{'lag_seconds'}, cellstr(sampleNames)]);
writetable(T_per, 'acf_per_sample_T-Windows_Ave_Example_1.csv');

% Mean ± SEM
T_mean = table(lag_seconds_vec, ACF_mean, ACF_sem, ...
    'VariableNames', {'lag_seconds','acf_mean','acf_sem'});
writetable(T_mean, 'acf_mean_T-Windows_Ave_Example_1.csv');

% ---------- Quick plot (optional) ----------
figure; hold on
plot(lag_seconds_vec, ACF, 'LineWidth', 0.6);
plot(lag_seconds_vec, ACF_mean, 'k', 'LineWidth', 2);
xlabel('Time lag (s)'); ylabel('Correlation coefficient');
title('Temporal ACF (per sample and mean)');
ylim([-0.2 1.05]); grid on;
legend([cellstr(sampleNames), {'mean'}], 'Location','bestoutside');

