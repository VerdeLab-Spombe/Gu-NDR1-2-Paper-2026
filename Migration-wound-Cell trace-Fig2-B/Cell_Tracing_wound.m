%% Cell migration tracing with fixed axis edges (ASCII-only)
% Fixed edges:
%   X max  = +500 um
%   Y min  = -150 um
%   Y max  = +200 um
% X min is auto-chosen to include your data (and 0), so the origin is visible.

clear; close all; clc

% ----- CONFIG -----
xlsxFile   = "Scr-shRNA-XY Coordinate-A1-1.xlsx";  % your file name
px2um      = 1.3043;                          % micrometers per pixel
figureID   = 1;

% Fixed axis edges (in micrometers)
xMaxFixed  = 200;   % X positive bound
yMinFixed  = -200;  % Y negative bound
yMaxFixed  = 200;   % Y positive bound
xMinFixed  = -200;

% ----- LOAD TRAJECTORIES (in pixels, relative to per-track origin) -----
[cellPos_px, origin_px, info] = loadXYTrajectories(xlsxFile);

% ----- CONVERT TO MICROMETERS -----
cellPos_um = cell(size(cellPos_px));
for ii = 1:numel(cellPos_px)
    P = cellPos_px{ii};
    if ~isempty(P)
        cellPos_um{ii} = P * px2um;   % scale relative coords
    else
        cellPos_um{ii} = [];
    end
end
origin_um = cellfun(@(o) o * px2um, origin_px, 'UniformOutput', false);

% ----- PLOT (units = micrometers) with fixed axis edges -----
plotTrajectoriesUM_FixedAxes(cellPos_um, figureID, xMaxFixed, yMinFixed, yMaxFixed, xMinFixed);

% ----- METRICS (units = micrometers, on relative traces) -----
metrics = computeMigrationMetrics(cellPos_um);
disp(metrics)
% writetable(metrics, 'cell_migration_metrics_um.csv');  % optional export


% ==================== FUNCTIONS ====================

function [cellPos, origin_px, info] = loadXYTrajectories(xlsxFile)
    % Reads X/Y columns from an Excel file and returns:
    %   cellPos: cell array of Nx2 doubles (relative traces in pixels)
    %   origin_px: cell array of 1x2 origins (first valid absolute point, pixels)
    % Detects X/Y by headers; falls back to 3-col [X,Y,extra] or 2-col [X,Y] blocks.

    T = readtable(xlsxFile, 'VariableNamingRule','preserve');
    ncol = size(T,2);
    vnames = T.Properties.VariableNames;

    % Try header-based X/Y pairs
    xyPairs = detectXYPairsByHeader(vnames);

    % Fallback 1: 3-col blocks after first column: [X, Y, extra]
    if isempty(xyPairs)
        if ncol >= 4
            cellNum = floor((ncol - 1)/3);
            tmp = zeros(cellNum,2);
            for k = 1:cellNum
                xidx = k*3;
                yidx = k*3 + 1;
                if yidx <= ncol
                    tmp(k,:) = [xidx, yidx];
                end
            end
            xyPairs = tmp(all(tmp>0,2),:);
        end
    end

    % Fallback 2: 2-col blocks after first column: [X, Y]
    if isempty(xyPairs)
        if ncol >= 3
            cellNum = floor((ncol - 1)/2);
            tmp = zeros(cellNum,2);
            for k = 1:cellNum
                xidx = 1 + (k-1)*2 + 1;  % skip first column (e.g., time)
                yidx = xidx + 1;
                if yidx <= ncol
                    tmp(k,:) = [xidx, yidx];
                end
            end
            xyPairs = tmp(all(tmp>0,2),:);
        end
    end

    cellPos   = {};
    origin_px = {};
    keep = false(size(xyPairs,1),1);

    for i = 1:size(xyPairs,1)
        cols = xyPairs(i,:);
        if any(cols < 1) || any(cols > ncol), continue, end

        Praw = table2array(T(:, cols));   % Nx2 absolute positions in pixels
        if isempty(Praw), continue, end

        % Remove rows with NaN in either X or Y
        Praw(any(isnan(Praw),2),:) = [];
        if isempty(Praw), continue, end

        % Origin = first valid absolute point (pixels)
        orig = Praw(1,:);
        Ptr  = Praw - orig;               % relative trace in pixels

        origin_px{end+1} = orig;          %#ok<AGROW>
        cellPos{end+1}   = Ptr;           %#ok<AGROW>
        keep(i) = true;
    end

    info = struct();
    info.VariableNames   = vnames;
    info.DetectedXYPairs = [];
    if any(keep)
        info.DetectedXYPairs = xyPairs(keep,:);
    end
end

function xyPairs = detectXYPairsByHeader(vnames)
    % Pair columns that look like base_X and base_Y (case-insensitive).
    % Also supports "X1","Y1" or "<base>_x<num>"/"<base>_y<num>".

    n = numel(vnames);
    cleaned = cellfun(@(s) lower(strtrim(regexprep(s, '\s+', ' '))), vnames, 'UniformOutput', false);

    bases = cell(n,1);
    isX = false(n,1);
    isY = false(n,1);

    for i = 1:n
        s = cleaned{i};

        % Ignore common non-position columns
        if contains(s,'time') || contains(s,'frame') || contains(s,'id')
            bases{i} = '';
            continue
        end

        % base + _x or _y (or space/dash)
        tok = regexp(s, '^(.*?)[ _\-]*([xy])$', 'tokens', 'once');
        if ~isempty(tok)
            bases{i} = strtrim(tok{1});
            if strcmp(tok{2},'x'), isX(i)=true; else, isY(i)=true; end
            continue
        end

        % xNN / yNN
        tok = regexp(s, '^([xy])(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            bases{i} = tok{2};
            if strcmp(tok{1},'x'), isX(i)=true; else, isY(i)=true; end
            continue
        end

        % base + _xNN / _yNN
        tok = regexp(s, '^(.*?)[ _\-]*([xy])(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            bases{i} = sprintf('%s_%s', strtrim(tok{1}), tok{3});
            if strcmp(tok{2},'x'), isX(i)=true; else, isY(i)=true; end
            continue
        end

        bases{i} = '';
    end

    xyPairs = [];
    if all(cellfun(@isempty, bases)), return, end

    [ubases, ~, ib] = unique(bases);
    for k = 1:numel(ubases)
        b = ubases{k};
        if isempty(b), continue, end
        idxs = find(ib == k);
        xidx = idxs(isX(idxs));
        yidx = idxs(isY(idxs));
        if numel(xidx) == 1 && numel(yidx) == 1
            xyPairs(end+1,:) = [xidx, yidx]; %#ok<AGROW>
        end
    end

    if ~isempty(xyPairs)
        xyPairs = sortrows(xyPairs,1);
    end
end

function plotTrajectoriesUM_FixedAxes(cellPos_um, figureID, xMaxFixed, yMinFixed, yMaxFixed)
    % Plot trajectories in micrometers; ticks with numeric labels.
    % Fixed edges: x upper bound, y lower & upper bounds. x lower bound auto.

    figure(figureID); clf
    set(gcf, 'Position', [520, 260, 560, 560]);
    ax = gca; hold(ax,'on'); box(ax,'on'); axis(ax,'equal');

    % Try to place axes at origin (R2020b+). Fallback draws zero lines.
    ok = true;
    try
        set(ax, 'XAxisLocation','origin', 'YAxisLocation','origin');
    catch
        ok = false;
    end
    if ~ok
        xline(0, ':');
        yline(0, ':');
    end

    % Gather data to choose xMin automatically (ensure origin visible)
    nonEmpty = ~cellfun(@isempty, cellPos_um);
    if any(nonEmpty)
        allXY = cell2mat(cellPos_um(nonEmpty)');
        xr = range(allXY(:,1)); yr = range(allXY(:,2));
        pad = 0.05 * max(1, max([xr, yr]));
        xMinAuto = min(allXY(:,1)) - pad;
        % Ensure origin is visible on the left if data are all positive
        xMin = min(xMinAuto, 0);
    else
        xMin = -500;   % fallback
    end

    % Apply fixed edges
    xlim([xMin, xMaxFixed]);
    ylim([yMinFixed, yMaxFixed]);

    % Warn if data exceed fixed limits (visible in Command Window)
    if any(nonEmpty)
        if any(allXY(:,1) > xMaxFixed)
            fprintf('Warning: some X data exceed xMaxFixed = %.3f um.\n', xMaxFixed);
        end
        if any(allXY(:,2) < yMinFixed) || any(allXY(:,2) > yMaxFixed)
            fprintf('Warning: some Y data exceed fixed Y limits [%.3f, %.3f] um.\n', yMinFixed, yMaxFixed);
        end
    end

    % Ticks with numeric labels
    targetTickCountX = 7;
    targetTickCountY = 7;
    [xt, xfmt] = niceTicks1D(xlim(ax), targetTickCountX);
    [yt, yfmt] = niceTicks1D(ylim(ax), targetTickCountY);
    set(ax, 'XTick', xt, 'YTick', yt);
    xtickformat(xfmt);
    ytickformat(yfmt);

    % Plot each track and its end point
    for ii = 1:numel(cellPos_um)
        P = cellPos_um{ii};
        if isempty(P), continue, end
        if size(P,1) >= 2
            plot(P(:,1), P(:,2), 'k-');
            plot(P(end,1), P(end,2), 'k.', 'MarkerSize', 18);
        else
            plot(P(1,1), P(1,2), 'k.', 'MarkerSize', 18);
        end
    end

    title('Cell Migration Traces (\mum)');
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
end

function [ticks, fmt] = niceTicks1D(lims, targetCount)
    % Compute "nice" tick positions and a matching numeric format.
    lo = min(lims); hi = max(lims);
    if ~isfinite(lo) || ~isfinite(hi) || lo == hi
        lo = lo - 1; hi = hi + 1;
    end
    rng = hi - lo;
    rough = rng / max(1, targetCount);
    mag = 10.^floor(log10(rough));
    bases = [1, 2, 5, 10];
    step = bases(find(bases.*mag >= rough, 1, 'first'));
    if isempty(step)
        step = 10*mag;
    else
        step = step * mag;
    end

    startTick = ceil(lo/step)*step;
    endTick   = floor(hi/step)*step;
    if startTick > endTick
        startTick = lo; endTick = hi;
    end
    ticks = startTick:step:endTick;

    dec = max(0, -floor(log10(step)));
    dec = min(dec, 3);
    fmt = sprintf('%%.%df', dec);
end

function M = computeMigrationMetrics(cellPos)
    % Metrics in same units as input (here, micrometers, relative traces).
    % N: time points; NetDisp: net displacement; PathLen: total path length;
    % Meander: NetDisp/PathLen; EndX/EndY: final coordinates.

    n = numel(cellPos);
    N = zeros(n,1);
    NetDisp = zeros(n,1);
    PathLen = zeros(n,1);
    Meander = zeros(n,1);
    EndX = zeros(n,1);
    EndY = zeros(n,1);

    for i = 1:n
        P = cellPos{i};
        if isempty(P)
            N(i) = 0; NetDisp(i) = NaN; PathLen(i) = NaN; Meander(i) = NaN; EndX(i) = NaN; EndY(i) = NaN;
            continue
        end
        N(i) = size(P,1);
        EndX(i) = P(end,1);
        EndY(i) = P(end,2);
        NetDisp(i) = hypot(EndX(i), EndY(i));
        if size(P,1) >= 2
            d = sqrt(sum(diff(P).^2, 2));
            PathLen(i) = sum(d);
            if PathLen(i) > 0
                Meander(i) = NetDisp(i) / PathLen(i);
            else
                Meander(i) = NaN;
            end
        else
            PathLen(i) = 0;
            Meander(i) = NaN;
        end
    end

    M = table((1:n)', N, NetDisp, PathLen, Meander, EndX, EndY, ...
        'VariableNames', {'CellID','N','NetDisp_um','PathLen_um','Meander','EndX_um','EndY_um'});
end
