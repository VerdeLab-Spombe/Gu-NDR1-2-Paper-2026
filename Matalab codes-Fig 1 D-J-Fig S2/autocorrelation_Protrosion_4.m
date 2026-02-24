%% Temporal ACF & dominant periodicity (CSV I/O)
% rows 2:26 = samples, cols 2:100 = time points (dt = 10 s)
clear; clc;

% -------- user params --------
infile   = 'Pro_Ave_example_1.csv';   % input CSV
dt_sec   = 10;      % sampling interval
lag_sec  = 200;     % total lag to each side (±200 s)
% Peak-picking settings (tune if needed)
minProm  = 0.10;    % minimum prominence to accept a peak (0..1)
minDistS = 2;       % minimum distance between peaks (in ACF samples)
% -----------------------------

% Read CSV and extract numeric block
M = readmatrix(infile);
%X = M(2:23, 2:101);     % rows = samples, cols = time
X = M; % for protrusions matrics
[numSamples, numTime] = size(X);     

% Clean & zero-mean each sample (row)
for i = 1:numSamples
    xi = X(i, :).';
    if any(isnan(xi))
        xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
    end
    xi = detrend(xi, 'constant');    % remove mean
    X(i, :) = xi.';
end

% Lags (±200 s at 10 s step => ±20 samples)
maxlag = round(lag_sec / dt_sec);     % 20
lag_samples_vec = (-maxlag:maxlag).';
lag_seconds_vec = lag_samples_vec * dt_sec;
nLags = numel(lag_seconds_vec);

% ACF per sample
ACF = zeros(nLags, numSamples);
for i = 1:numSamples
    [acfi, ~] = xcorr(X(i, :).', maxlag, 'coeff');
    ACF(:, i) = acfi(:);
end

% Mean ± SEM across samples
ACF_mean = mean(ACF, 2, 'omitnan');
ACF_sem  = std(ACF, 0, 2, 'omitnan') ./ sqrt(numSamples);

% -------- Dominant periodicity (per sample) --------
posMask       = lag_samples_vec > 0;            % positive lags only
lagsPosSamp   = lag_samples_vec(posMask);
lagsPosSec    = lag_seconds_vec(posMask);

peakLagSamples = nan(numSamples,1);
peakLagSeconds = nan(numSamples,1);
peakValue      = nan(numSamples,1);
peakProm       = nan(numSamples,1);

for i = 1:numSamples
    y = ACF(posMask, i);
    [idx, pk, prom] = dominantPeak(y, minProm, minDistS);
    if ~isnan(idx)
        peakLagSamples(i) = lagsPosSamp(idx);
        peakLagSeconds(i) = lagsPosSec(idx);
        peakValue(i)      = pk;
        peakProm(i)       = prom;
    end
end

% Dominant periodicity for the space-mean ACF (optional but useful)
yMean = ACF_mean(posMask);
[idxM, pkM, promM] = dominantPeak(yMean, minProm, minDistS);
meanACF_period_sec   = nan;
meanACF_period_samps = nan;
if ~isnan(idxM)
    meanACF_period_sec   = lagsPosSec(idxM);
    meanACF_period_samps = lagsPosSamp(idxM);
end

% Summary stats across samples (ignore NaNs)
valid = ~isnan(peakLagSeconds);
n_valid = sum(valid);
meanPeriod = mean(peakLagSeconds(valid));
medianPeriod = median(peakLagSeconds(valid));
stdPeriod = std(peakLagSeconds(valid));
semPeriod = stdPeriod / sqrt(max(n_valid,1));

% --------------- Write CSV outputs ---------------
% Lags
T_lags = table(lag_seconds_vec, lag_samples_vec, ...
    'VariableNames', {'lag_seconds','lag_samples'});
writetable(T_lags, 'Pro_Ave_example_1_lags.csv');

% Per-sample ACF
sampleNames = "Sample_" + string(1:numSamples);
T_per = array2table([lag_seconds_vec, ACF], ...
    'VariableNames', [{'lag_seconds'}, cellstr(sampleNames)]);
writetable(T_per, 'Pro_Ave_example_1_acf_per_sample.csv');

% Mean ± SEM ACF
T_mean = table(lag_seconds_vec, ACF_mean, ACF_sem, ...
    'VariableNames', {'lag_seconds','acf_mean','acf_sem'});
writetable(T_mean, 'Pro_Ave_example_1_acf_mean.csv');

% Per-sample dominant periodicity
T_perPeriod = table( (1:numSamples).', peakLagSeconds, peakLagSamples, ...
                     peakValue, peakProm, ...
    'VariableNames', {'sample_id','period_seconds','period_samples','peak_value','peak_prominence'});
writetable(T_perPeriod, 'Pro_Ave_example_1_periodicity_per_sample.csv');

% Summary of periodicity
T_sum = table( n_valid, meanPeriod, medianPeriod, stdPeriod, semPeriod, ...
               meanACF_period_sec, meanACF_period_samps, dt_sec, lag_sec, minProm, minDistS, ...
    'VariableNames', {'n_valid','mean_period_s','median_period_s','std_period_s','sem_period_s', ...
                      'meanACF_peak_period_s','meanACF_peak_period_samples', ...
                      'dt_seconds','max_lag_seconds','min_prominence','min_peak_distance_samples'});
writetable(T_sum, 'Pro_Ave_example_1_periodicity_summary.csv');

% --------------- (Optional) quick plot ---------------
figure; hold on
plot(lag_seconds_vec, ACF, 'LineWidth', 0.6);
plot(lag_seconds_vec, ACF_mean, 'k', 'LineWidth', 2);
xlabel('Time lag (s)'); ylabel('Correlation coefficient');
title('Temporal ACF (per sample + mean)');
ylim([-0.2 1.05]); grid on;
legend([cellstr(sampleNames), {'mean'}], 'Location','bestoutside');

%% -------- helper: dominant peak finder (excludes zero lag) --------
function [idx, pk, prom] = dominantPeak(y, minProm, minDist)
% Return location (index into y), value, and prominence of the
% dominant (highest) positive-lag peak. If no valid peak, returns NaNs.
    idx = NaN; pk = NaN; prom = NaN;
    % Prefer findpeaks if available
    if exist('findpeaks','file') == 2
        [pks, locs, ~, proms] = findpeaks(y, ...
            'MinPeakProminence', minProm, 'MinPeakDistance', minDist);
        if ~isempty(pks)
            [pk, k] = max(pks);
            idx  = locs(k);
            prom = proms(k);
        end
        return;
    end
    % Fallback: naive local-maximum detection
    n = numel(y);
    if n < 3, return; end
    isPeak = [false; y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end); false];
    candIdx = find(isPeak);
    if isempty(candIdx), return; end
    [pk, k] = max(y(candIdx));
    idx = candIdx(k);
    prom = NaN; % not computed in fallback
end
